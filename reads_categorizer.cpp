#include <iostream>
#include <algorithm>
#include <numeric>

#include "utils.h"
#include "sam_utils.h"
#include "libs/cptl_stl.h"

std::mutex mtx;
std::string workspace;
config_t config;

std::mutex* mtx_contig;
std::vector<std::vector<std::string> > mate_seqs;

std::vector<uint32_t> depths;
uint64_t qual_counts[256];
std::vector<std::vector<uint32_t> > dist_between_end_and_rnd;

samFile* get_writer(std::string dir, std::string name, bam_hdr_t* header) {
    samFile* writer = sam_open((dir + "/" + name).c_str(), "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + (workspace + name);
    }
    return writer;
}

// computes the differences for the left and the right half of the read
std::pair<int, int> compute_left_and_right_differences(bam1_t* r) {
    int border = r->core.l_qseq/2;
    int left_Ms = 0;

    // compute how many Ms (either matches or mismatches are in the first half of the read)
    // in other words, this is readlen/2 - number of insertions in the first half
    int qpos = 0;
    uint32_t* cigar = bam_get_cigar(r);
    for (int i = 0; i < r->core.n_cigar; i++) {
        char op_char = bam_cigar_opchr(cigar[i]);
        int op_len = bam_cigar_oplen(cigar[i]);
        if (op_char == 'M') {
            left_Ms += std::min(border, qpos+op_len) - qpos;
            qpos += op_len;
        } else if (op_char == 'I') {
            qpos += op_len;
        }

        if (qpos >= border) break;
    }

    // we computed left_Ms because the MD tag does not include insertions
    // i.e. if I have CIGAR: 50M50I50M and MD: 74A25, the mismatch is not in position 75 in the query,
    // but it is in position 125
    // here we compute the number of deleted bases and mismatches in each half of the read based on the MD tag
    std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
    int m = 0;
    int left_diffs = 0, right_diffs = 0;
    bool del_mode = false;
    for (int i = 0; i < md_tag.length(); i++) {
        char c = md_tag[i];
        if (c >= '0' && c <= '9') {
            m = m*10 + c-'0';
            del_mode = false;
        } else if (c >= 'A' && c <= 'T') {
            left_Ms -= m;
            m = 0;
            if (left_Ms > 0) left_diffs++;
            else right_diffs++;
            if (!del_mode) left_Ms--;
        } else if (c == '^') {
            del_mode = true;
        }
    }

    // We also need to count the number of inserted bases, since the MD tag does not report this information
    qpos = 0;
    for (int i = 0; i < r->core.n_cigar; i++) {
        char op_char = bam_cigar_opchr(cigar[i]);
        int op_len = bam_cigar_oplen(cigar[i]);
        if (op_char == 'M') {
            qpos += op_len;
        } else if (op_char == 'I') {
            if (qpos < border && qpos + op_len > border) {  // this is the case there an insertion is partially in the left
                                                            // half and partially in the right half
                left_diffs += border - qpos;
                right_diffs += qpos + op_len - border;
            } else if (qpos < border) {
                left_diffs += op_len;
            } else if (qpos >= border) {
                right_diffs += op_len;
            }
            qpos += op_len;
        }
    }

    return {left_diffs, right_diffs};
}

const int MIN_HSR_DIFFS = 4;
// we are being very conservative, we require at least 4 differences in one half and perfect match in the other half
bool is_left_hsr(bam1_t* r) {
	auto p = compute_left_and_right_differences(r);
	return !is_left_clipped(r, 0) && !is_right_clipped(r, 0) && p.first >= MIN_HSR_DIFFS && p.second == 0;
}
bool is_right_hsr(bam1_t* r) {
	auto p = compute_left_and_right_differences(r);
	return !is_left_clipped(r, 0) && !is_right_clipped(r, 0) && p.first == 0 && p.second >= MIN_HSR_DIFFS;
}

int find_left_M_seq(bam1_t* r) { // assumes read is not clipped
	std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
	int m1 = 0;
	for (int i = 0; i < md_tag.length() && md_tag[i] >= '0' && md_tag[i] <= '9'; i++) {
		m1 = m1*10 + (md_tag[i]-'0');
	}
	int m2 = bam_cigar_oplen(bam_get_cigar(r)[0]); // guaranteed to be M
	return std::min(m1, m2);
}
int find_right_M_seq(bam1_t* r) {
	std::string md_tag = bam_aux2Z(bam_aux_get(r, "MD"));
	int m1 = 0;
	int d = 1;
	for (int i = md_tag.length()-1; i >= 0 && md_tag[i] >= '0' && md_tag[i] <= '9'; i--) {
		m1 += d*(md_tag[i]-'0');
		d *= 10;
	}
	int m2 = bam_cigar_oplen(bam_get_cigar(r)[r->core.n_cigar-1]); // guaranteed to be M
	return std::min(m1, m2);
}

void update_cigar(bam1_t* read, uint32_t* cigar, int n_cigar) {
	int l_aux = bam_get_l_aux(read);
	int l_data = read->core.l_qname + 4*n_cigar + (read->core.l_qseq+1)/2
			+ read->core.l_qseq + l_aux;
	uint32_t m_data = l_data;
	kroundup32(m_data);
	uint8_t* data = new uint8_t[m_data];
	memset(data, 0, m_data);

	uint8_t* mov_data = data;
	memcpy(mov_data, (uint8_t*) bam_get_qname(read), read->core.l_qname);
	mov_data += read->core.l_qname;
	memcpy(mov_data, cigar, 4*n_cigar);
	mov_data += 4*n_cigar;

	memcpy(mov_data, (uint8_t*) bam_get_seq(read), (read->core.l_qseq+1)/2);
	mov_data += (read->core.l_qseq+1)/2;
	memcpy(mov_data, (uint8_t*) bam_get_qual(read), read->core.l_qseq);
	mov_data += read->core.l_qseq;
	memcpy(mov_data, (uint8_t*) bam_get_aux(read), l_aux);

	read->l_data = l_data;
	read->m_data = m_data;
	read->data = data;
	read->core.n_cigar = n_cigar;
}

void categorize(int id, int contig_id, std::string contig_name, std::string bam_fname, std::string reference_fname,
		std::vector<int> rnd_positions) {
    mtx.lock();
    std::cout << "Categorizing " << contig_name << std::endl;
    mtx.unlock();

    std::sort(rnd_positions.begin(), rnd_positions.end());

    open_samFile_t* bam_file = open_samFile(bam_fname);
    if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
        throw "Failed to read reference " + reference_fname;
    }

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    if (iter == NULL) { // no reads
    	close_samFile(bam_file);
    	return;
    }

    bam1_t* read = bam_init1();

    samFile* clip_writer = NULL;
    samFile* rdc_writer = NULL;
	samFile* ldc_writer = NULL;

    int curr_pos = 0;
	std::vector<uint32_t> rnd_positions_depth(rnd_positions.size());

	uint64_t contig_qual_counts[256];
	std::fill(contig_qual_counts, contig_qual_counts+256, uint64_t(0));

	std::vector<std::vector<uint32_t> > rnd_positions_dist_between_end_and_rnd(rnd_positions.size());

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (!is_primary(read)) continue;

		int64_t mq = get_mq(read);
		// note that unmapped reads do not have a MQ tag, so we specifically have to let them in
        if (is_dc_pair(read) && (read->core.qual >= config.min_stable_mapq || mq >= config.min_stable_mapq || is_unmapped(read))) {
        	if (read->core.qual >= mq && !is_unmapped(read)) { // stable end
        		if (bam_is_rev(read) && !is_right_clipped(read, config.min_clip_len)) {
        			if (!ldc_writer) ldc_writer = get_writer(workspace + "/L", std::to_string(contig_id) + ".noremap.bam", bam_file->header);
        			int ok = sam_write1(ldc_writer, bam_file->header, read);
					if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        		} else if (!bam_is_rev(read) && !is_left_clipped(read, config.min_clip_len)) {
        			if (!rdc_writer) rdc_writer = get_writer(workspace + "/R", std::to_string(contig_id) + ".noremap.bam", bam_file->header);
        			int ok = sam_write1(rdc_writer, bam_file->header, read);
					if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        		}
        	}
        	if (read->core.qual <= mq || is_unmapped(read)) { // save read seq for remapping
        		std::string qname = bam_get_qname(read), read_seq = get_sequence(read, true), qual_ascii = get_qual_ascii(read, true);
				if (is_samechr(read)) {
					if (read->core.flag & BAM_FREAD1) qname += "_1";
					else qname += "_2";
				}
				mtx_contig[read->core.mtid].lock();
				mate_seqs[read->core.mtid].push_back(qname + " " + read_seq + " " + qual_ascii);
				mtx_contig[read->core.mtid].unlock();
        	}
        }

        if (is_unmapped(read)) continue;

		while (curr_pos < rnd_positions.size() && read->core.pos > rnd_positions[curr_pos]) curr_pos++;

		// sample depth
		bool sampled = false;
		hts_pos_t read_endpos = bam_endpos(read);
		for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < read_endpos; i++) {
			if (read->core.pos <= rnd_positions[i] && rnd_positions[i] <= read_endpos) {
				sampled = true;
				rnd_positions_depth[i]++;
			}
		}

		// sample pairs crossing
		hts_pos_t pair_startpos = read->core.pos + config.read_len/2, pair_endpos = get_mate_endpos(read) - config.read_len/2;
		if (read->core.isize > 0 && read->core.isize <= config.max_is && is_samechr(read) && !is_samestr(read)
			&& !is_left_clipped(read, config.min_clip_len) && !is_right_clipped(read, config.min_clip_len)) {
			for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < pair_endpos; i++) {
				 if (pair_startpos <= rnd_positions[i] && rnd_positions[i] <= pair_endpos) {
					rnd_positions_dist_between_end_and_rnd[i].push_back(pair_endpos-rnd_positions[i]);
				}
			}
		}

		// TODO: HSRs do not seem to have a big impact on insertions
//		if (is_right_hsr(read)) {
//			int m = find_left_M_seq(read);
//			uint32_t o1 = bam_cigar_gen(m, BAM_CMATCH), o2 = bam_cigar_gen(read->core.l_qseq-m, BAM_CSOFT_CLIP);
//			uint32_t c[2] = { o1, o2 };
//			update_cigar(read, c, 2);
//			std::string m_str = std::to_string(m);
//			bam_aux_update_str(read, "MD", m_str.length(), m_str.c_str());
//		} else if (is_left_hsr(read)) {
//			int m = find_right_M_seq(read);
//			hts_pos_t end = bam_endpos(read);
//			read->core.pos = end - m;
//			uint32_t o1 = bam_cigar_gen(read->core.l_qseq-m, BAM_CSOFT_CLIP), o2 = bam_cigar_gen(m, BAM_CMATCH);
//			uint32_t c[2] = { o1, o2 };
//			update_cigar(read, c, 2);
//			std::string m_str = std::to_string(m);
//			bam_aux_update_str(read, "MD", m_str.length(), m_str.c_str());
//		}

		if (is_left_clipped(read, config.min_clip_len) || is_right_clipped(read, config.min_clip_len)) {
			if (!clip_writer) clip_writer = get_writer(workspace + "/clipped/", std::to_string(contig_id) + ".bam", bam_file->header);

			int ok = sam_write1(clip_writer, bam_file->header, read);
			if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
		} else if (sampled) {
			contig_qual_counts[int(avg_qual(read)+0.5)]++;
		}
    }

    if (clip_writer) sam_close(clip_writer);
    if (rdc_writer) sam_close(rdc_writer);
    if (ldc_writer) sam_close(ldc_writer);

    mtx.lock();
	for (uint32_t d : rnd_positions_depth) {
		if (d > 0) depths.push_back(d);
	}
	for (int i = 0; i < 256; i++) qual_counts[i] += contig_qual_counts[i];
	for (auto& v : rnd_positions_dist_between_end_and_rnd) {
		if (v.empty()) continue;
		dist_between_end_and_rnd.push_back(v);
	}
	mtx.unlock();

	bam_destroy1(read);
	hts_itr_destroy(iter);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    std::string workdir = std::string(argv[2]);
    workspace = workdir + "/workspace/";
    std::string reference_fname = argv[3];

    contig_map_t contig_map;
    contig_map.parse(workdir);
    config.parse(workdir + "/config.txt");

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());
	if (hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str())) != 0) {
		throw "Failed to read reference " + reference_fname;
	}
	bam_hdr_t* header = bam_file->header;

    mtx_contig = new std::mutex[header->n_targets];
    mate_seqs.resize(header->n_targets);

    // read random positions
	std::string contig_name;
	std::ifstream rnd_pos_fin(workdir + "/random_pos.txt");
	int pos;
	std::unordered_map<std::string, std::vector<int> > rnd_pos_map;
	while (rnd_pos_fin >> contig_name >> pos) {
		rnd_pos_map[contig_name].push_back(pos);
	}

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool.push(categorize, contig_id, contig_name, bam_fname, reference_fname, rnd_pos_map[contig_name]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        }
    }

    for (int i = 0; i < header->n_targets; i++) {
    	if (mate_seqs[i].empty()) continue;
		std::string fname = std::to_string(i) + ".txt";
		std::ofstream mate_seqs_fout(workspace + "/mateseqs/" + fname);
		for (std::string& mate_seq : mate_seqs[i]) {
			mate_seqs_fout << mate_seq << std::endl;
		}
		mate_seqs_fout.close();
	}

    // get 1-percentile of base qualities
    uint64_t tot_quals = std::accumulate(qual_counts, qual_counts+256, uint64_t(0));
    uint64_t _1_perc_count = tot_quals/100;
    int min_avg_base_qual = 0;
    for (int i = 0; i < 256; i++) {
    	if (_1_perc_count <= qual_counts[i]) {
    		min_avg_base_qual = i;
    		break;
    	}
    	_1_perc_count -= qual_counts[i];
    }

    // for each insertion size, compute the estimate the minimum number of discordant pairs
    // in theory, it would be the number of discordant pairs *per breakpoint*, so we would need to require twice as many from the insertion
    // in practice, we do not because insertions often show a reduced number of pairs, perhaps due to overrepresentation of low complexity
    // subsequences in inserted sequences
    std::vector<std::vector<uint32_t>> pairs_crossing_dists(config.max_is+1);
	for (int i = 0; i < dist_between_end_and_rnd.size(); i++) {
		std::vector<uint32_t> dist(config.max_is+1);
		for (uint32_t val : dist_between_end_and_rnd[i]) {
			if (val <= config.max_is) dist[val]++;
		}
		pairs_crossing_dists[0].push_back(dist[0]);
		for (int j = 1; j <= config.max_is; j++) {
			dist[j] += dist[j-1];
			pairs_crossing_dists[j].push_back(dist[j]);
		}
	}

	std::ofstream mdpbs_fout(workdir + "/min_disc_pairs_by_size.txt");
	for (int i = 0; i <= config.max_is; i++) {
		std::sort(pairs_crossing_dists[i].begin(), pairs_crossing_dists[i].end());
		mdpbs_fout << i << " " << pairs_crossing_dists[i][pairs_crossing_dists[i].size()/200] << std::endl;
	}
	mdpbs_fout.close();

    std::ofstream stats_out(workdir + "/stats.txt");
	std::sort(depths.begin(), depths.end());
	stats_out << "min_depth " << depths[depths.size()/100] << std::endl;
	stats_out << "median_depth " << depths[depths.size()/2] << std::endl;
	stats_out << "max_depth " << depths[depths.size()-depths.size()/100] << std::endl;
    stats_out << "min_avg_base_qual " << min_avg_base_qual << std::endl;
	stats_out.close();
}
