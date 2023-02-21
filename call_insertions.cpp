#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>

#include "libs/ssw_cpp.h"
#include "vcf_utils.h"
#include "sam_utils.h"
#include "utils.h"
#include "libs/cptl_stl.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

std::mutex mtx;
config_t config;
std::string workspace;

chr_seqs_map_t contigs;

std::vector<insertion_t*> insertions;

int MAX_BP_DIST = 10;

struct clip_consensus_t {
    hts_pos_t pos;
    std::string seq;
    int fwd_clipped, rev_clipped;

    clip_consensus_t(hts_pos_t pos, std::string& seq, int fwd_clipped, int rev_clipped) :
    	pos(pos), seq(seq), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped) {}
};

insertion_t* get_insertion(std::string& contig_name, clip_consensus_t& rc_consensus, clip_consensus_t& lc_consensus,
                           StripedSmithWaterman::Aligner& aligner) {
    suffix_prefix_aln_t spa = aln_suffix_prefix(rc_consensus.seq, lc_consensus.seq, 1, -4, config.max_seq_error);
    if (spa.overlap < config.min_clip_len || is_homopolymer(lc_consensus.seq.c_str(), spa.overlap)) return NULL;

    int mm = spa.mismatches;
    int aln_len = spa.overlap;
    double mm_rate = double(mm)/aln_len;
    if (mm_rate > config.max_seq_error) return NULL;

    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment aln_rh, aln_lh;

    std::string consensus_junction_seq = std::string(config.clip_penalty, 'N') + rc_consensus.seq +
    		lc_consensus.seq.substr(spa.overlap) + std::string(config.clip_penalty, 'N');
    hts_pos_t offset = rc_consensus.pos - consensus_junction_seq.length();
    if (offset < 0) offset = 0;
    aligner.Align(consensus_junction_seq.c_str(), contigs.get_seq(contig_name) + offset,
                  consensus_junction_seq.length(), filter, &aln_lh, 0);
    if (aln_lh.query_begin > 0 || aln_lh.query_end-aln_lh.query_begin+1 < config.clip_penalty+config.min_clip_len) return NULL;

    hts_pos_t left_bp = offset + aln_lh.ref_end;
    std::string left_anchor = std::to_string(offset + aln_lh.ref_begin) + "-" + std::to_string(offset + aln_lh.ref_end);

    std::string remaining_consensus = consensus_junction_seq.substr(aln_lh.query_end+1);
    offset = rc_consensus.pos;
    hts_pos_t len = consensus_junction_seq.length();
    if (offset + len > contigs.get_len(contig_name)) len = contigs.get_len(contig_name) - offset;
    aligner.Align(remaining_consensus.c_str(), contigs.get_seq(contig_name) + offset,
                  len, filter, &aln_rh, 0);
    if (aln_rh.query_end < remaining_consensus.length()-1 || aln_rh.query_end-aln_rh.query_begin+1 < config.clip_penalty+config.min_clip_len) return NULL;

    hts_pos_t right_bp = offset + aln_rh.ref_begin-1;
    std::string right_anchor = std::to_string(offset + aln_rh.ref_begin) + "-" + std::to_string(offset + aln_rh.ref_end);

    hts_pos_t ins_seq_len = consensus_junction_seq.length() - (aln_lh.query_end+1) - (remaining_consensus.length()-aln_rh.query_begin);
    std::string ins_seq = consensus_junction_seq.substr(aln_lh.query_end+1, ins_seq_len);
    if (ins_seq_len == 0) return NULL;

    insertion_t* insertion = new insertion_t(contig_name, left_bp, right_bp, 0, 0,
			rc_consensus.fwd_clipped, rc_consensus.rev_clipped, lc_consensus.fwd_clipped, lc_consensus.rev_clipped, spa.overlap, ins_seq);
    insertion->left_anchor = left_anchor, insertion->right_anchor = right_anchor;

    return insertion;
}

void call_insertions(int id, int contig_id, std::string contig_name) {
    mtx.lock();
    std::cout << "Calling insertions for " << contig_name << std::endl;
    mtx.unlock();

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, true);

    std::ifstream clip_fin(workspace + "/clip_consensuses/" + std::to_string(contig_id) + ".txt");
    std::vector<clip_consensus_t> rc_consensuses, lc_consensuses;
    std::string chr, dir, seq;
    hts_pos_t start, end, breakpoint;
    int fwd_clipped, rev_clipped;
    while (clip_fin >> chr >> start >> end >> breakpoint >> dir >> seq >> fwd_clipped >> rev_clipped) {
        if (dir == "R") {
            rc_consensuses.push_back(clip_consensus_t(breakpoint, seq, fwd_clipped, rev_clipped));
        } else {
            lc_consensuses.push_back(clip_consensus_t(breakpoint, seq, fwd_clipped, rev_clipped));
        }
    }

    for (int i = 0; i < rc_consensuses.size(); i++) {
        clip_consensus_t& rc_consensus = rc_consensuses[i];
        for (int j = 0; j < lc_consensuses.size(); j++) {
            clip_consensus_t& lc_consensus = lc_consensuses[j];
            int mh_len = rc_consensus.pos - lc_consensus.pos;
            if (lc_consensus.pos-rc_consensus.pos <= MAX_BP_DIST && mh_len <= config.max_mh_len) {
                insertion_t* insertion = get_insertion(contig_name, rc_consensus, lc_consensus, aligner);
                if (insertion != NULL && mh_len < (int) insertion->ins_seq.length()) {
                    mtx.lock();
                    insertions.push_back(insertion);
                    mtx.unlock();
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {

    std::string workdir = argv[1];
    workspace = workdir + "/workspace/";
    std::string reference_fname = argv[2];
    std::string sample_name = argv[3];

    std::string full_cmd_fname = workdir + "/full_cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    contigs.read_fasta_into_map(reference_fname);

    contig_map_t contig_map;
    contig_map.parse(workdir);
    config.parse(workdir + "/config.txt");

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = thread_pool.push(call_insertions, contig_id, contig_name);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    std::string out_vcf_fname = workdir + "/small_ins.vcf.gz";
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
	}

	bcf_hdr_t* out_vcf_header = generate_vcf_header(contigs, sample_name, config.min_insertion_size, full_cmd_str, config.version);
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	std::sort(insertions.begin(), insertions.end(), [&out_vcf_header](insertion_t* i1, insertion_t* i2) {
		int contig_id1 = bcf_hdr_name2id(out_vcf_header, i1->chr.c_str());
		int contig_id2 = bcf_hdr_name2id(out_vcf_header, i2->chr.c_str());
		// negative because we want descending order
		int sc_score1 = -(i1->rc_reads()*i1->lc_reads()), sc_score2 = -(i2->rc_reads()*i2->lc_reads());
		int overlap1 = -i1->overlap, overlap2 = -i2->overlap;
		return std::tie(contig_id1, i1->start, i1->end, i1->ins_seq, sc_score2, overlap1) <
			   std::tie(contig_id2, i2->start, i2->end, i2->ins_seq, sc_score1, overlap2);
	});

    int int_id = 0;
    bcf1_t* bcf_entry = bcf_init();
    std::unordered_set<std::string> used_keys;
    for (insertion_t* insertion : insertions) {
    	std::string key = unique_key(insertion);
    	if (used_keys.count(key)) continue;
    	used_keys.insert(key);

    	std::string id = "S_INS_" + std::to_string(int_id);
		insertion_to_bcf_entry(insertion, out_vcf_header, bcf_entry, id, contigs);

		bcf_update_info_int32(out_vcf_header, bcf_entry, "OVERLAP", &insertion->overlap, 1);
		bcf_update_info_string(out_vcf_header, bcf_entry, "ALGORITHM", "consensus_overlap");
		bcf_update_info_string(out_vcf_header, bcf_entry, "LEFT_ANCHOR", insertion->left_anchor.c_str());
		bcf_update_info_string(out_vcf_header, bcf_entry, "RIGHT_ANCHOR", insertion->right_anchor.c_str());

        if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
        	throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
        }
        int_id++;
    }

    bcf_close(out_vcf_file);
}
