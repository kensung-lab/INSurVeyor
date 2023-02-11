#include <iostream>
#include <queue>
#include <unordered_map>
#include <numeric>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "cluster.h"
#include "libs/cptl_stl.h"
#include "utils.h"
#include "vcf_utils.h"

std::string workdir;
std::mutex mtx;

stats_t stats;
config_t config;
contig_map_t contig_map;

int MAX_READ_IS;

const int MAX_BUFFER_SIZE = 100;

std::queue<open_samFile_t*> bam_pool;

open_samFile_t* get_bam_reader(std::string& bam_fname, std::string& reference_fname) {
    if (!bam_pool.empty()) {
        open_samFile_t* o = bam_pool.front();
        bam_pool.pop();
        return o;
    }

    open_samFile_t* o = open_samFile(bam_fname.c_str());
    hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
    return o;
}

void release_bam_reader(open_samFile_t* reader) {
    bam_pool.push(reader);
}


void find_spanning(int id, insertion_t* insertion, std::string bam_fname, std::string reference_fname) {
    mtx.lock();
    open_samFile_t* bam_file = get_bam_reader(bam_fname, reference_fname);
    mtx.unlock();

    char r_region[1000], l_region[1000];
    std::stringstream l_region_ss, r_region_ss;
	l_region_ss << insertion->chr << ":" << std::max(hts_pos_t(1), insertion->start-config.max_is) << "-" << insertion->start+config.max_is;
	r_region_ss << insertion->chr << ":" << std::max(hts_pos_t(1), insertion->end-config.max_is) << "-" << insertion->end+config.max_is;
	strcpy(l_region, l_region_ss.str().c_str());
	strcpy(r_region, r_region_ss.str().c_str());

	char* regions[] = {l_region, r_region};

	std::unordered_map<std::string, int> qname_is_concordant;

	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read) || !(read->core.flag & BAM_FPROPER_PAIR)) continue;
        if (is_left_clipped(read, config.min_clip_len) || is_right_clipped(read, config.min_clip_len)) continue;
        if (avg_qual(read) < stats.get_min_avg_base_qual()) continue;

        if (read->core.isize > 0) {
			int pair_start = read->core.pos+config.read_len/2, pair_end = get_mate_endpos(read)-config.read_len/2;
			int b = 0;
			if (pair_start <= insertion->start && insertion->start <= pair_end) b |= 1;
			if (pair_start <= insertion->end && insertion->end <= pair_end) b |= 2;
			qname_is_concordant[bam_get_qname(read)] = b;
        } else if (read->core.isize < 0) {
        	int b = qname_is_concordant[bam_get_qname(read)];
        	if (b & 1) insertion->l_conc_pairs++;
        	if (b & 2) insertion->r_conc_pairs++;
        }
    }

    mtx.lock();
    release_bam_reader(bam_file);
    mtx.unlock();
};

void compute_coverage(int id, insertion_t* insertion, std::string bam_fname, std::string reference_fname) {
	mtx.lock();
	open_samFile_t* bam_file = get_bam_reader(bam_fname, reference_fname);
	mtx.unlock();

	char c;
	std::stringstream ss_la(insertion->left_anchor), ss_ra(insertion->right_anchor);
	hts_pos_t left_flanking_start, left_flanking_end, right_flanking_start, right_flanking_end;
	ss_la >> left_flanking_start >> c >> left_flanking_end;
	ss_ra >> right_flanking_start >> c >> right_flanking_end;
	hts_pos_t left_flanking_len = left_flanking_end - left_flanking_start;
	hts_pos_t right_flanking_len = right_flanking_end - right_flanking_start;

	std::stringstream l_region_ss, r_region_ss;
	l_region_ss << insertion->chr << ":" << left_flanking_start << "-" << left_flanking_end;
	r_region_ss << insertion->chr << ":" << right_flanking_start << "-" << right_flanking_end;
	char l_region[1000], r_region[1000];
	strcpy(l_region, l_region_ss.str().c_str());
	strcpy(r_region, r_region_ss.str().c_str());

	char* regions[] = {l_region, r_region};

	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
	bam1_t* read = bam_init1();
	std::vector<int64_t> flanking_left_cov(left_flanking_len), flanking_right_cov(right_flanking_len);
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		hts_pos_t rs = read->core.pos, re = bam_endpos(read);

		hts_pos_t b = std::max(hts_pos_t(0), rs-left_flanking_start), e = std::min(left_flanking_len, re-left_flanking_start);
		for (hts_pos_t i = b; i < e; i++) {
			flanking_left_cov[i]++;
		}

		b = std::max(hts_pos_t(0), rs-right_flanking_start), e = std::min(right_flanking_len, re-right_flanking_start);
		for (hts_pos_t i = b; i < e; i++) {
			flanking_right_cov[i]++;
		}
	}

	std::sort(flanking_left_cov.begin(), flanking_left_cov.end(), std::greater<int64_t>());
	int end = std::find(flanking_left_cov.begin(), flanking_left_cov.end(), 0)-flanking_left_cov.begin();
	insertion->median_lf_cov = flanking_left_cov[end/2];

	std::sort(flanking_right_cov.begin(), flanking_right_cov.end(), std::greater<int64_t>());
	end = std::find(flanking_right_cov.begin(), flanking_right_cov.end(), 0)-flanking_right_cov.begin();
	insertion->median_rf_cov = flanking_right_cov[end/2];

	mtx.lock();
    release_bam_reader(bam_file);
    mtx.unlock();
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    std::string reference_fname = argv[3];

    config.parse(workdir + "/config.txt");
    contig_map.parse(workdir);
    stats.parse_stats(workdir + "/stats.txt", config.per_contig_stats);

    std::vector<std::pair<insertion_t*, bcf1_t*> > small_insertions, assembled_insertions, transurveyor_insertions;

    bcf1_t* bcf_entry = bcf_init();

    std::string in_small_ins_fname = workdir + "/small_ins.vcf.gz";
    htsFile* in_small_ins_file = bcf_open(in_small_ins_fname.c_str(), "r");
    bcf_hdr_t* small_ins_hdr = bcf_hdr_read(in_small_ins_file);
    while (bcf_read(in_small_ins_file, small_ins_hdr, bcf_entry) == 0) {
    	insertion_t* ins = new insertion_t(bcf_seqname_safe(small_ins_hdr, bcf_entry), bcf_entry->pos,
    			get_sv_end(bcf_entry, small_ins_hdr), 0, 0, 0, 0, 0, "");
    	ins->left_anchor = get_left_anchor(bcf_entry, small_ins_hdr);
    	ins->right_anchor = get_right_anchor(bcf_entry, small_ins_hdr);
    	small_insertions.push_back({ins, bcf_dup(bcf_entry)});
    }

    std::string in_assembled_ins_fname = workdir + "/assembled_ins.vcf.gz";
    htsFile* in_assembled_ins_file = bcf_open(in_assembled_ins_fname.c_str(), "r");
	bcf_hdr_t* assembled_ins_hdr = bcf_hdr_read(in_assembled_ins_file);
	while (bcf_read(in_assembled_ins_file, assembled_ins_hdr, bcf_entry) == 0) {
		insertion_t* ins = new insertion_t(bcf_seqname_safe(assembled_ins_hdr, bcf_entry), bcf_entry->pos,
				get_sv_end(bcf_entry, assembled_ins_hdr), 0, 0, 0, 0, 0, "");
		ins->left_anchor = get_left_anchor(bcf_entry, assembled_ins_hdr);
		ins->right_anchor = get_right_anchor(bcf_entry, assembled_ins_hdr);
		assembled_insertions.push_back({ins, bcf_dup(bcf_entry)});
	}

	std::string in_transurveyor_ins_fname = workdir + "/transurveyor_ins.vcf.gz";
	htsFile* in_transurveyor_ins_file = bcf_open(in_transurveyor_ins_fname.c_str(), "r");
	bcf_hdr_t* transurveyor_ins_hdr = bcf_hdr_read(in_transurveyor_ins_file);
	while (bcf_read(in_transurveyor_ins_file, transurveyor_ins_hdr, bcf_entry) == 0) {
		insertion_t* ins = new insertion_t(bcf_seqname_safe(transurveyor_ins_hdr, bcf_entry), bcf_entry->pos,
				get_sv_end(bcf_entry, transurveyor_ins_hdr), 0, 0, 0, 0, 0, "");
		ins->left_anchor = get_left_anchor(bcf_entry, transurveyor_ins_hdr);
		ins->right_anchor = get_right_anchor(bcf_entry, transurveyor_ins_hdr);
		int* discordant = NULL, * split_reads = NULL;
		int n = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "DISCORDANT", &discordant, &n);
		n = 0;
		bcf_get_info_int32(transurveyor_ins_hdr, bcf_entry, "SPLIT_READS", &split_reads, &n);
		if (discordant[0]+split_reads[0] > 1 && discordant[1]+split_reads[1] > 1) {
			transurveyor_insertions.push_back({ins, bcf_dup(bcf_entry)});
		}
	}

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (auto& p : small_insertions) {
    	std::future<void> future = thread_pool.push(compute_coverage, p.first, bam_fname, reference_fname);
    	futures.push_back(std::move(future));
    }
    for (auto& p : assembled_insertions) {
		std::future<void> future1 = thread_pool.push(compute_coverage, p.first, bam_fname, reference_fname);
		futures.push_back(std::move(future1));
		std::future<void> future2 = thread_pool.push(find_spanning, p.first, bam_fname, reference_fname);
		futures.push_back(std::move(future2));
	}
    for (auto& p : transurveyor_insertions) {
		std::future<void> future1 = thread_pool.push(compute_coverage, p.first, bam_fname, reference_fname);
		futures.push_back(std::move(future1));
		std::future<void> future2 = thread_pool.push(find_spanning, p.first, bam_fname, reference_fname);
		futures.push_back(std::move(future2));
	}
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    std::string out_small_ins_fname = workdir + "/small_ins.annotated.vcf.gz";
	htsFile* out_small_ins_file = bcf_open(out_small_ins_fname.c_str(), "wz");
	if (bcf_hdr_write(out_small_ins_file, small_ins_hdr) != 0) {
		throw std::runtime_error("Failed to write header to " + out_small_ins_fname);
	}
	for (auto& p : small_insertions) {
		int int2_conv[] = {p.first->median_lf_cov, p.first->median_rf_cov};
		bcf_update_info_int32(small_ins_hdr, p.second, "STABLE_DEPTHS", int2_conv,  2);
		if (bcf_write(out_small_ins_file, small_ins_hdr, p.second) != 0) {
			throw std::runtime_error("Failed to write VCF record to " + out_small_ins_fname);
		}
	}

	std::string out_assembled_ins_fname = workdir + "/assembled_ins.annotated.vcf.gz";
	htsFile* out_assembled_ins_file = bcf_open(out_assembled_ins_fname.c_str(), "wz");
	if (bcf_hdr_write(out_assembled_ins_file, assembled_ins_hdr) != 0) {
		throw std::runtime_error("Failed to write header to " + out_assembled_ins_fname);
	}
	for (auto& p : assembled_insertions) {
		int int2_conv[] = {p.first->median_lf_cov, p.first->median_rf_cov};
		bcf_update_info_int32(assembled_ins_hdr, p.second, "STABLE_DEPTHS", int2_conv,  2);
		int2_conv[0] = p.first->l_conc_pairs; int2_conv[1] = p.first->r_conc_pairs;
		bcf_update_info_int32(assembled_ins_hdr, p.second, "SPANNING_READS", int2_conv,  2);
		if (bcf_write(out_assembled_ins_file, assembled_ins_hdr, p.second) != 0) {
			throw std::runtime_error("Failed to write VCF record to " + out_assembled_ins_fname);
		}
	}

	std::string out_transurveyor_ins_fname = workdir + "/transurveyor_ins.annotated.vcf.gz";
	htsFile* out_transurveyor_ins_file = bcf_open(out_transurveyor_ins_fname.c_str(), "wz");
	if (bcf_hdr_write(out_transurveyor_ins_file, transurveyor_ins_hdr) != 0) {
		throw std::runtime_error("Failed to write header to " + out_transurveyor_ins_fname);
	}
	for (auto& p : transurveyor_insertions) {
		int int2_conv[] = {p.first->median_lf_cov, p.first->median_rf_cov};
		bcf_update_info_int32(transurveyor_ins_hdr, p.second, "STABLE_DEPTHS", int2_conv,  2);
		int2_conv[0] = p.first->l_conc_pairs; int2_conv[1] = p.first->r_conc_pairs;
		bcf_update_info_int32(transurveyor_ins_hdr, p.second, "SPANNING_READS", int2_conv,  2);
		if (bcf_write(out_transurveyor_ins_file, transurveyor_ins_hdr, p.second) != 0) {
			throw std::runtime_error("Failed to write VCF record to " + out_transurveyor_ins_fname);
		}
	}

	bcf_close(in_small_ins_file);
	bcf_close(out_small_ins_file);
	bcf_close(in_assembled_ins_file);
	bcf_close(out_assembled_ins_file);
	bcf_close(in_transurveyor_ins_file);
	bcf_close(out_transurveyor_ins_file);

    while (!bam_pool.empty()) {
        close_samFile(bam_pool.front());
        bam_pool.pop();
    }
}
