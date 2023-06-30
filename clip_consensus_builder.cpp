#include <iostream>
#include <algorithm>
#include <cmath>

#include "sam_utils.h"
#include "utils.h"
#include "libs/cptl_stl.h"

std::mutex mtx;
config_t config;
std::string workdir, workspace;

struct consensus_t {
    bool left_clipped;
    int contig_id;
    hts_pos_t start, breakpoint, end;
    std::string consensus;
    int fwd_clipped, rev_clipped;

    consensus_t(bool left_clipped, int contig_id, hts_pos_t start, hts_pos_t breakpoint, hts_pos_t end,
                const std::string& consensus, int fwd_clipped, int rev_clipped)
                : left_clipped(left_clipped), contig_id(contig_id), start(start), breakpoint(breakpoint), end(end),
                consensus(consensus), fwd_clipped(fwd_clipped), rev_clipped(rev_clipped) {}
};

struct del_ins_t {
    int del, ins;

    del_ins_t() : del(0), ins(0) {};
};
del_ins_t get_dels_ins_in_first_n_chars(bam_redux_t* r, int n) {
    del_ins_t del_ins;
    int offset = 0;
    for (uint32_t c : r->cigar) {
        int len = bam_cigar_oplen(c);
        char op = bam_cigar_opchr(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = bam_cigar_type(c) & 2 || op == 'S';
        if (consumes_ref && offset + len > n) {
            len = n-offset;
        }

        if (op == 'D') {
            del_ins.del += len;
        } else if (op == 'I') {
            del_ins.ins += len;
        }

        if (consumes_ref) {
            offset += len;
            if (offset == n) break;
        }
    }
    return del_ins;
}

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};
bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.freq != bs2.freq) return bs1.freq < bs2.freq;
    return bs1.qual < bs2.qual;
}

hts_pos_t get_start_offset(bam_redux_t* r1, bam_redux_t* r2) {
    hts_pos_t offset = r2->unclipped_start() - r1->unclipped_start();
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    del_ins_t del_ins = get_dels_ins_in_first_n_chars(r1, offset);
    return offset + del_ins.ins - del_ins.del;
}

std::string build_full_consensus_seq(std::vector<bam_redux_t*>& clipped) {
    const int MAX_CONSENSUS_LEN = 100000;
    char consensus[MAX_CONSENSUS_LEN];

    std::vector<hts_pos_t> read_start_offsets, read_end_offsets;
    hts_pos_t consensus_len = 0;
    for (bam_redux_t* r : clipped) {
        hts_pos_t start_offset = get_start_offset(clipped[0], r);
        read_start_offsets.push_back(start_offset);

        hts_pos_t end_offset = start_offset + r->seq_len() - 1;
        read_end_offsets.push_back(end_offset);
        consensus_len = std::max(consensus_len, end_offset+1);
    }

    int s = 0;
    for (int i = 0; i < consensus_len; i++) {
        while (s < clipped.size() && read_end_offsets[s] < i) s++;

        base_score_t base_scores[4] = { base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T') };
        for (int j = s; j < clipped.size() && read_start_offsets[j] <= i; j++) {
            if (read_end_offsets[j] < i) continue;

            char nucl = get_base(clipped[j]->seq.data(), i - read_start_offsets[j]);
            uint8_t qual = clipped[j]->qual[i - read_start_offsets[j]];
            if (nucl == 'A') {
                base_scores[0].freq++;
                base_scores[0].qual += qual;
            } else if (nucl == 'C') {
                base_scores[1].freq++;
                base_scores[1].qual += qual;
            } else if (nucl == 'G') {
                base_scores[2].freq++;
                base_scores[2].qual += qual;
            } else if (nucl == 'T') {
                base_scores[3].freq++;
                base_scores[3].qual += qual;
            }
        }

        base_score_t best_base_score = max(base_scores[0], base_scores[1], base_scores[2], base_scores[3]);

        consensus[i] = best_base_score.base;
    }
    consensus[consensus_len] = '\0';

    return std::string(consensus);
}

consensus_t* build_full_consensus(int contig_id, std::vector<bam_redux_t*>& clipped, bool left_clipped) {
    if (clipped.size() <= 2) return NULL;

    std::sort(clipped.begin(), clipped.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
        return r1->unclipped_start() < r2->unclipped_start();
    });

    std::string consensus_seq = build_full_consensus_seq(clipped);
    if (consensus_seq == "") return NULL;

    std::vector<bam_redux_t*> accepted_reads;
    for (bam_redux_t* r : clipped) {
        int mm = 0;
        int mm_clip = 0; // mismatches in the clipped portion only

        // filter reads with too many differences from the consensus_seq
        hts_pos_t offset = get_start_offset(clipped[0], r);
        hts_pos_t clip_start = left_clipped ? 0 : r->seq_len() - r->right_clip_size;
        hts_pos_t clip_end = left_clipped ? r->left_clip_size : r->seq_len();
        for (int i = 0; i < r->seq_len(); i++) {
            if (i+offset >= consensus_seq.length()) {
                std::cerr << "WARNING: consensus_seq out of boundary." << std::endl;
            }
            if (consensus_seq[i + offset] != get_base(r->seq.data(), i)) {
                mm++;
                if (i >= clip_start && i < clip_end) mm_clip++;
            }
        }
        if (mm <= std::ceil(config.max_seq_error * consensus_seq.length()) &&
        mm_clip <= (clip_end-clip_start)/2) {
            /* this is meant to fix a corner case where the clip is very short, and completely different from
             * the consensus_seq. However, the rest of the reads is the same as the consensus_seq, therefore the mismatches
             * accumulated in the clip are not enough to discard the read, despite the read not belonging to the cluster.
             * We are being very permissive (i.e. we allow up to 50% mismatches in the clip) because
             * 1. Illumina is known to accumulate sequencing errors at the 3' tail, which is the clipped portion in about
             * 50% of the cases
             * 2. Especially for short clips, if we used config.max_seq_error as for the rest of the consensus_seq, 2-3 mismatches
             * would already discard them
             */
            accepted_reads.push_back(r);
        }
    }

    clipped.swap(accepted_reads);
    if (clipped.size() <= 2) return NULL;

    hts_pos_t start = clipped[0]->unclipped_start(), end = 0;
    hts_pos_t breakpoint = left_clipped ? INT32_MAX : 0; // the current HTS_POS_MAX does not compile on some compilers
    int fwd_clipped = 0, rev_clipped = 0;
    for (bam_redux_t* r : clipped) {
        breakpoint = left_clipped ? std::min(breakpoint, r->start) : std::max(breakpoint, r->end);
        end = std::max(end, r->unclipped_end());
        if (r->is_rev()) rev_clipped++;
        else fwd_clipped++;
    }

    consensus_seq = build_full_consensus_seq(clipped);

    // these bps have support from only one or two reads, so they are prone to errors
    hts_pos_t remove_from_start = get_start_offset(clipped[0], clipped[2]);
    sort(clipped.begin(), clipped.end(), [](bam_redux_t* r1, bam_redux_t* r2) {
        return r1->unclipped_end() > r2->unclipped_end();
    });
    hts_pos_t remove_from_end = clipped[0]->unclipped_end() - clipped[2]->unclipped_end();

    consensus_seq = consensus_seq.substr(remove_from_start, consensus_seq.length() - remove_from_start - remove_from_end);
    consensus_t* consensus = new consensus_t(left_clipped, contig_id, start+remove_from_start, breakpoint,
                                             end-remove_from_end, consensus_seq, fwd_clipped, rev_clipped);
    return consensus;
}

void build_consensuses(int id, int contig_id, std::string contig_name, hts_pos_t contig_len) {

    std::string clip_fname = workspace + "/clipped/" + std::to_string(contig_id) + ".bam";
    if (!file_exists(clip_fname)) return;

    mtx.lock();
    std::cout << "Building consensuses for " << contig_name << std::endl;
    mtx.unlock();

    open_samFile_t* bam_file = open_samFile(clip_fname, true);
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::vector<bam_redux_t*> lc_reads, rc_reads;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_left_clipped(read, config.min_clip_len) && !is_right_clipped(read, config.min_clip_len)) {
        	bam_redux_t* bam_redux = new bam_redux_t(read);
        	if (bam_redux->unclipped_end() < contig_len) {
				lc_reads.push_back(bam_redux);
        	} else {
        		delete bam_redux;
        	}
        }
        if (is_right_clipped(read, config.min_clip_len) && !is_left_clipped(read, config.min_clip_len)) {
        	bam_redux_t* bam_redux = new bam_redux_t(read);
        	if (bam_redux->unclipped_start() >= 0) {
				rc_reads.push_back(bam_redux);
        	} else {
        		delete bam_redux;
        	}
        }
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
    close_samFile(bam_file);

    std::vector<consensus_t*> full_consensuses;
    std::vector<bam_redux_t*> curr_candidate_cluster;

    auto lc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->start-r2->start) <= config.max_clipped_pos_dist;};
    if (!lc_reads.empty()) {
        for (bam_redux_t* lc_read : lc_reads) {
            if (!curr_candidate_cluster.empty() &&
            !lc_same_cluster(curr_candidate_cluster[0], lc_read)) { // candidate cluster complete
                consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, true);
                curr_candidate_cluster.clear();
                if (full_consensus != NULL) {
                    full_consensuses.push_back(full_consensus);
                }
            }
            curr_candidate_cluster.push_back(lc_read);
        }
        // process last cluster
        consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, true);
        curr_candidate_cluster.clear();
        if (full_consensus != NULL) {
            full_consensuses.push_back(full_consensus);
        }
    }
    auto rc_same_cluster = [](bam_redux_t* r1, bam_redux_t* r2) {return abs(r1->end-r2->end) <= config.max_clipped_pos_dist;};
    sort(rc_reads.begin(), rc_reads.end(), [](bam_redux_t* r1, bam_redux_t* r2) { return r1->end < r2->end; });
    if (!rc_reads.empty()) {
        for (bam_redux_t* rc_read : rc_reads) {
            if (!curr_candidate_cluster.empty() &&
            !rc_same_cluster(curr_candidate_cluster[0], rc_read)) { // candidate cluster complete
                consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, false);
                if (full_consensus != NULL) {
                    full_consensuses.push_back(full_consensus);
                }
                curr_candidate_cluster.clear();
            }
            curr_candidate_cluster.push_back(rc_read);
        }
        // process last cluster
        consensus_t* full_consensus = build_full_consensus(contig_id, curr_candidate_cluster, false);
        curr_candidate_cluster.clear();
        if (full_consensus != NULL) {
            full_consensuses.push_back(full_consensus);
        }
    }

    for (bam_redux_t* r : lc_reads) delete r;
    for (bam_redux_t* r : rc_reads) delete r;

    if (full_consensuses.empty()) return;

    std::ofstream clip_fout(workspace + "/clip_consensuses/" + std::to_string(contig_id) + ".txt");
    for (consensus_t* consensus : full_consensuses) {
        clip_fout << contig_name << " " << consensus->start << " " << consensus->end << " " << consensus->breakpoint << " ";
        clip_fout << (consensus->left_clipped ? "L " : "R ") << consensus->consensus << " " << consensus->fwd_clipped << " " << consensus->rev_clipped << std::endl;
    }
    clip_fout.close();
}

int main(int argc, char* argv[]) {
    workdir = argv[1];
    workspace = workdir + "/workspace";

    std::string reference_fname = argv[2];

    contig_map_t contig_map;
    contig_map.parse(workdir);
    config.parse(workdir + "/config.txt");

    chr_seqs_map_t contigs;
    contigs.read_fasta_into_map(reference_fname);

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        hts_pos_t contig_len = contigs.get_len(contig_name);
        std::future<void> future = thread_pool.push(build_consensuses, contig_id, contig_name, contig_len);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();
}
