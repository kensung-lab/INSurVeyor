#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_

#include <mutex>

#include "utils.h"
#include "dc_remapper.h"
#include "libs/ssw.h"

extern stats_t stats;
extern config_t config;
extern std::ofstream graph_writer;
extern std::ofstream assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer,
assembly_failed_bad_anchors_writer, assembly_failed_mh_too_long, assembly_failed_lt50bp, assembly_succeeded;

extern const int TOO_MANY_READS;

std::mutex graph_mtx, failed_assembly_mtx;

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

struct path_permission_t {
	bool can_start_path, can_end_path;
};
struct seq_w_pp_t {
	std::string seq;
	path_permission_t clip_pair;

	seq_w_pp_t() : seq(), clip_pair() {}
	seq_w_pp_t(std::string& seq, bool can_start_path, bool can_end_path) : seq(seq) {
		clip_pair.can_start_path = can_start_path;
		clip_pair.can_end_path = can_end_path;
	}
};

std::vector<std::string> assemble_reads(std::vector<seq_w_pp_t>& left_stable_read_seqs, std::vector<seq_w_pp_t>& unstable_read_seqs,
		std::vector<seq_w_pp_t>& right_stable_read_seqs, StripedSmithWaterman::Aligner& aligner, config_t config, std::stringstream& ss_graph);

std::vector<int> find_rev_topological_order(int n, std::vector<int>& out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

void build_aln_guided_graph(std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& alns, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev) {
	std::sort(alns.begin(), alns.end(),
			[](const std::pair<std::string, StripedSmithWaterman::Alignment>& aln1, const std::pair<std::string, StripedSmithWaterman::Alignment>& aln2) {
		return aln1.second.ref_begin < aln2.second.ref_begin;
	});

	for (int i = 0; i < alns.size(); i++) {
		for (int j = i+1; j < alns.size() && alns[i].second.ref_end-alns[j].second.ref_begin >= config.min_clip_len; j++) {
			suffix_prefix_aln_t spa = aln_suffix_prefix(alns[i].first, alns[j].first, 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa.overlap) {
				out_edges[i]++;
				l_adj[i].push_back({j, spa.score, spa.overlap});
				l_adj_rev[j].push_back({i, spa.score, spa.overlap});
			}
		}
	}

	// two major differences with the regular assembly:
	// 1 - by how the graph is defined, no cycle is possible here
	// 2 - in regular assembly, we only report contigs made of at last 2 reads. Here we report even single reads
}

bool accept(const StripedSmithWaterman::Alignment& aln, std::string& log, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	// the 'N' we insterted for padding are counted as mismatches, so they must be removed from the count
	int lc_size = get_left_clip_size(aln), rc_size = get_right_clip_size(aln);
	if (lc_size > config.clip_penalty && rc_size > config.clip_penalty) return false; // do not accept if both sides clipped

	std::stringstream ss;
	bool lc_is_noise = false, rc_is_noise = false;
	if (!qual_ascii.empty()) {
		if (lc_size > config.clip_penalty && lc_size-config.clip_penalty < qual_ascii.length()/2) {
			std::string lc_qual = qual_ascii.substr(0, lc_size-config.clip_penalty);
			std::string not_lc_qual = qual_ascii.substr(lc_size-config.clip_penalty);
			lc_is_noise = avg_qual(lc_qual) < min_avg_base_qual && avg_qual(not_lc_qual) >= min_avg_base_qual;
			ss << aln.cigar_string << " " << lc_qual << " " << avg_qual(lc_qual) << " " << not_lc_qual << " " << avg_qual(not_lc_qual) << " " << lc_is_noise << std::endl;
		}
		if (rc_size > config.clip_penalty && rc_size-config.clip_penalty < qual_ascii.length()/2) {
			std::string not_rc_qual = qual_ascii.substr(0, qual_ascii.length()-(rc_size-config.clip_penalty));
			std::string rc_qual = qual_ascii.substr(qual_ascii.length()-(rc_size-config.clip_penalty));
			rc_is_noise = avg_qual(rc_qual) < min_avg_base_qual && avg_qual(not_rc_qual) >= min_avg_base_qual;
			ss << aln.cigar_string << " " << rc_qual << " " << avg_qual(rc_qual) << " " << not_rc_qual << " " << avg_qual(not_rc_qual) << " " << rc_is_noise << std::endl;
		}
	}
	log += ss.str();

	int left_padding_matched = std::max(0, config.clip_penalty - lc_size);
	int right_padding_matched = std::max(0, config.clip_penalty - rc_size);
	int mismatches = aln.mismatches - left_padding_matched - right_padding_matched;
	double mismatch_rate = double(mismatches)/(aln.query_end-aln.query_begin-left_padding_matched-right_padding_matched);
	return (lc_size <= config.clip_penalty || lc_is_noise) && (rc_size <= config.clip_penalty || rc_is_noise) && mismatch_rate <= max_seq_error;
};
bool accept(const StripedSmithWaterman::Alignment& aln, double max_seq_error = 1.0, std::string qual_ascii = "", int min_avg_base_qual = 0) {
	std::string temp;
	return accept(aln, temp, max_seq_error, qual_ascii, min_avg_base_qual);
}
void add_alignment(std::string& reference, std::string& query, std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& accepted_alns,
		std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> >& rejected_alns, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	std::string padded_query = std::string(config.clip_penalty, 'N') + query + std::string(config.clip_penalty, 'N');
	aligner.Align(padded_query.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
	aln.ref_begin += config.clip_penalty - get_left_clip_size(aln);
	aln.ref_end -= config.clip_penalty - get_right_clip_size(aln);
	if (accept(aln)) {
		accepted_alns.push_back({query, aln});
	} else {
		rejected_alns.push_back({query, aln});
	}
}

void correct_contig(std::string& contig, std::vector<std::string>& reads, StripedSmithWaterman::Aligner& harsh_aligner) {
	std::vector<int> As(contig.length()), Cs(contig.length()), Gs(contig.length()), Ts(contig.length());
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	for (std::string& read : reads) {
		std::string padded_read = std::string(config.clip_penalty, 'N') + read + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_read.c_str(), contig.c_str(), contig.length(), filter, &aln, 0);
		if (accept(aln)) {
			int left_padding_aligned = config.clip_penalty - get_left_clip_size(aln);
			for (int i = 0; i < read.length(); i++) {
				char c = read[i];
				if (c == 'A') As[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'C') Cs[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'G') Gs[i+aln.ref_begin+left_padding_aligned]++;
				else if (c == 'T') Ts[i+aln.ref_begin+left_padding_aligned]++;
			}
		}
	}

	for (int i = 0; i < contig.length(); i++) {
		int max_freq = max(As[i], Cs[i], Gs[i], Ts[i]);
		if (max_freq == 0) continue;
		if (max_freq == As[i]) contig[i] = 'A';
		else if (max_freq == Cs[i]) contig[i] = 'C';
		else if (max_freq == Gs[i]) contig[i] = 'G';
		else if (max_freq == Ts[i]) contig[i] = 'T';
	}
}

std::vector<std::string> generate_reference_guided_consensus(std::string reference, reads_cluster_t* r_cluster, reads_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<StripedSmithWaterman::Alignment>& consensus_contigs_alns, std::string& consensus_log) {

	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > accepted_alns, _rejected_alns_lf, _rejected_alns_is, _rejected_alns_rf;
	for (bam1_t* read : r_cluster->reads) {
		std::string read_seq = get_sequence(read);
		add_alignment(reference, read_seq, accepted_alns, _rejected_alns_lf, aligner);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		rc(mate_seq);
		add_alignment(reference, mate_seq, accepted_alns, _rejected_alns_is, aligner);
	}
	for (bam1_t* read : l_cluster->reads) {
		std::string read_seq = get_sequence(read);
		add_alignment(reference, read_seq, accepted_alns, _rejected_alns_rf, aligner);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		add_alignment(reference, mate_seq, accepted_alns, _rejected_alns_is, aligner);
	}
	if (r_cluster->clip_cluster) {
		add_alignment(reference, r_cluster->clip_cluster->full_seq, accepted_alns, _rejected_alns_lf, aligner);
	}
	if (l_cluster->clip_cluster) {
		add_alignment(reference, l_cluster->clip_cluster->full_seq, accepted_alns, _rejected_alns_rf, aligner);
	}

	int n = accepted_alns.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);
	build_aln_guided_graph(accepted_alns, out_edges, l_adj, l_adj_rev);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = accepted_alns[curr_vertex].first;
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(accepted_alns[curr_vertex].first);
		while (best_edges[curr_vertex].overlap) {
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += accepted_alns[curr_vertex].first.substr(overlap);
			used_reads.push_back(accepted_alns[curr_vertex].first);
		}
		used[curr_vertex] = true;

		std::string corrected_assembled_sequence = assembled_sequence;
		correct_contig(corrected_assembled_sequence, used_reads, harsh_aligner);
		assembled_sequences.push_back(corrected_assembled_sequence);
	}
	for (int i = 0; i < n; i++) {
		if (!used[i]) assembled_sequences.push_back(accepted_alns[i].first);
	}

	// retain assembled sequences that align without clipping and do not overlap a higher rated sequence
	std::vector<std::string> retained_assembled_sequences;
	for (std::string& assembled_sequence : assembled_sequences) {
		aligner.Align(assembled_sequence.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln)) continue;

		bool overlaps = false;
		for (StripedSmithWaterman::Alignment& existing_aln : consensus_contigs_alns) {
			if (std::max(existing_aln.ref_begin, aln.ref_begin) <= std::min(existing_aln.ref_end, aln.ref_end)) {
				overlaps = true;
				break;
			}
		}
		if (!overlaps) {
			consensus_contigs_alns.push_back(aln);
			retained_assembled_sequences.push_back(assembled_sequence);
		}
	}

	if (retained_assembled_sequences.empty()) return {};

	int l_bp_in_seq = r_cluster->end()-r_cluster->start(), r_bp_in_seq = reference.length()-(l_cluster->end()-l_cluster->start());

	std::stringstream ss_res2, ss;

	// try scaffolding using rejected reads
	std::stringstream ss_graph;
	std::vector<seq_w_pp_t> rejected_alns_lf, rejected_alns_is, rejected_alns_rf;
	for (auto& e : _rejected_alns_lf) rejected_alns_lf.push_back({e.first, true, true});
	for (auto& e : _rejected_alns_is) rejected_alns_is.push_back({e.first, true, true});
	for (auto& e : _rejected_alns_rf) rejected_alns_rf.push_back({e.first, true, true});
	std::vector<char> rejected_alns_lf_clipped(rejected_alns_lf.size(), 'N'), rejected_alns_is_clipped(rejected_alns_is.size(), 'N'),
			rejected_alns_rf_clipped(rejected_alns_rf.size(), 'N');
	std::vector<std::string> scaffolds = assemble_reads(rejected_alns_lf, rejected_alns_is, rejected_alns_rf,
			harsh_aligner, config, ss_graph);

	ss << "Assembled " << n << " reads into " << assembled_sequences.size() << " sequences." << std::endl;
	ss << reference.length() << " " << l_bp_in_seq << " " << r_bp_in_seq << std::endl;
	for (std::string a : assembled_sequences) {
		aligner.Align(a.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		ss << a.length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << " ";
	}
	ss << std::endl;
	for (int i = 0; i < retained_assembled_sequences.size(); i++) {
		StripedSmithWaterman::Alignment& aln = consensus_contigs_alns[i];
		ss << retained_assembled_sequences[i].length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << " ";
	}
	ss << std::endl << std::endl;

	ss << "RETAINED CONTIGS: " << retained_assembled_sequences.size() <<  std::endl;
	for (std::string& a : retained_assembled_sequences) ss << a << std::endl;
	ss << "SCAFFOLDS: " << scaffolds.size() << std::endl;
	for (std::string& a : scaffolds) ss << a << std::endl;
	ss << std::endl;

	std::vector<std::pair<std::string, StripedSmithWaterman::Alignment> > contigs_sorted_by_pos;
	for (int i = 0; i < consensus_contigs_alns.size(); i++) {
		contigs_sorted_by_pos.push_back({retained_assembled_sequences[i], consensus_contigs_alns[i]});
	}
	std::sort(contigs_sorted_by_pos.begin(), contigs_sorted_by_pos.end(),
			[](std::pair<std::string, StripedSmithWaterman::Alignment>& p1, std::pair<std::string, StripedSmithWaterman::Alignment>& p2) {
		return p1.second.ref_begin < p2.second.ref_begin;
	});

	std::vector<int> linked(contigs_sorted_by_pos.size()-1, -1);
	std::vector<std::pair<int, int> > link_overlaps(contigs_sorted_by_pos.size()-1);
	for (int i = 0; i < scaffolds.size(); i++) {
		std::string& scaffold = scaffolds[i];
		int best_link = -1, best_link_w = 0;
		std::pair<int, int> best_link_overlap;
		for (int j = 0; j < contigs_sorted_by_pos.size()-1; j++) {
			if (linked[j] >= 0) continue;

			suffix_prefix_aln_t spa1 = aln_suffix_prefix(contigs_sorted_by_pos[j].first, scaffold, 1, -4, config.max_seq_error, config.min_clip_len);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(scaffold, contigs_sorted_by_pos[j+1].first, 1, -4, config.max_seq_error, config.min_clip_len);
			if (spa1.overlap && spa2.overlap && best_link_w < spa1.score+spa2.score) {
				best_link = j, best_link_w = spa1.score+spa2.score;
				best_link_overlap = {spa1.overlap, spa2.overlap};
			}
		}

		if (best_link >= 0) {
			linked[best_link] = i;
			link_overlaps[best_link] = best_link_overlap;
		}
	}

	std::vector<std::string> scaffolded_sequences;
	std::string curr_seq = contigs_sorted_by_pos[0].first;
	for (int i = 1; i < contigs_sorted_by_pos.size(); i++) {
		if (linked[i-1] == -1) {
			scaffolded_sequences.push_back(curr_seq);
			curr_seq = contigs_sorted_by_pos[i].first;
		} else {
			std::string link = scaffolds[linked[i-1]];
			auto& lo = link_overlaps[i-1];
			link = link.substr(lo.first, link.length()-lo.first-lo.second);
			curr_seq += link + contigs_sorted_by_pos[i].first;
		}
	}
	scaffolded_sequences.push_back(curr_seq);

	std::vector<StripedSmithWaterman::Alignment> scaffolded_seqs_alns;
	ss << "SCAFFOLDED SEQUENCES: " << scaffolded_sequences.size() << std::endl;
	bool scaffolding_failed = false;
	for (std::string s : scaffolded_sequences) {
		aligner.Align(s.c_str(), reference.c_str(), reference.length(), filter, &aln, 0);
		if (!accept(aln)) {
			scaffolding_failed = true;
			break;
		}
		scaffolded_seqs_alns.push_back(aln);
		ss << s.length() << "," << aln.ref_begin << "-" << aln.ref_end << "," << accept(aln) << " ";
	} ss << std::endl << std::endl;

	if (!scaffolding_failed)
	for (int i = 0; i < scaffolded_seqs_alns.size(); i++) {
		for (int j = i+1; j < scaffolded_seqs_alns.size(); j++) {
			if (std::max(scaffolded_seqs_alns[i].ref_begin, scaffolded_seqs_alns[j].ref_begin) <= std::min(scaffolded_seqs_alns[i].ref_end, scaffolded_seqs_alns[j].ref_end)) {
				scaffolding_failed = true;
				break;
			}
		}
		if (scaffolding_failed) break;
	}

	ss << reference << std::endl;

	consensus_log = ss_res2.str() + ss.str();

	if (!scaffolding_failed) scaffolded_seqs_alns.swap(consensus_contigs_alns);
	return scaffolding_failed ? retained_assembled_sequences : scaffolded_sequences;
}

void build_graph(std::vector<std::string>& read_seqs, std::vector<int>& order, std::vector<int>& out_edges,
		std::vector<std::vector<edge_t> >& l_adj, std::vector<std::vector<edge_t> >& l_adj_rev,
		int max_mismatches, int min_overlap, std::stringstream& ss_graph) {

	int n = read_seqs.size();

	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			std::string& s1 = read_seqs[i];
			std::string& s2 = read_seqs[j];

			int max_overlap = std::min(s1.length(), s2.length())-1;
			suffix_prefix_aln_t spa1 = aln_suffix_prefix(s1, s2, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa1_homopolymer = is_homopolymer(s2.c_str(), spa1.overlap);
			suffix_prefix_aln_t spa2 = aln_suffix_prefix(s2, s1, 1, -4, 1.0, min_overlap, max_overlap, max_mismatches);
			bool spa2_homopolymer = is_homopolymer(s1.c_str(), spa2.overlap);
			if (spa1.overlap && spa2.overlap) {
				if (spa1.score >= spa2.score && order[i] <= order[j]) {
					ss_graph << i << " -> " << j << " " << spa1.score << " " << spa1.overlap << " " << spa1_homopolymer << std::endl;
					if (!spa1_homopolymer) {
						out_edges[i]++;
						l_adj[i].push_back({j, spa1.score, spa1.overlap});
						l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
					}
				} else if (spa1.score < spa2.score && order[j] <= order[i]) {
					ss_graph << j << " -> " << i << " " << spa2.score << " " << spa2.overlap << " " << spa2_homopolymer << std::endl;
					if (!spa2_homopolymer) {
						out_edges[j]++;
						l_adj[j].push_back({i, spa2.score, spa2.overlap});
						l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
					}
				}
			} else if (spa1.overlap && order[i] <= order[j]) {
				ss_graph << i << " -> " << j << " " << spa1.score << " " << spa1.overlap << " " << spa1_homopolymer << std::endl;
				if (!spa1_homopolymer) {
					out_edges[i]++;
					l_adj[i].push_back({j, spa1.score, spa1.overlap});
					l_adj_rev[j].push_back({i, spa1.score, spa1.overlap});
				}
			} else if (spa2.overlap && order[j] <= order[i]) {
				ss_graph << j << " -> " << i << " " << spa2.score << " " << spa2.overlap << " " << spa2_homopolymer << std::endl;
				if (!spa2_homopolymer) {
					out_edges[j]++;
					l_adj[j].push_back({i, spa2.score, spa2.overlap});
					l_adj_rev[i].push_back({j, spa2.score, spa2.overlap});
				}
			}
		}
	} ss_graph << std::endl;
}

std::vector<std::string> assemble_reads(std::vector<seq_w_pp_t>& left_stable_read_seqs, std::vector<seq_w_pp_t>& unstable_read_seqs,
		std::vector<seq_w_pp_t>& right_stable_read_seqs, StripedSmithWaterman::Aligner& harsh_aligner, config_t config, std::stringstream& ss_graph) {

	std::vector<std::string> read_seqs;
	std::vector<path_permission_t> path_permissions;
	std::vector<int> order;
	for (seq_w_pp_t& s : left_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(1);
	}
	for (seq_w_pp_t& s : unstable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(2);
	}
	for (seq_w_pp_t& s : right_stable_read_seqs) {
		read_seqs.push_back(s.seq);
		path_permissions.push_back(s.clip_pair);
		order.push_back(3);
	}

	int n = read_seqs.size();
	std::vector<int> out_edges(n);
	std::vector<std::vector<edge_t> > l_adj(n), l_adj_rev(n);

	for (int i = 0; i < n; i++) {
		ss_graph << i << " " << read_seqs[i] << " " << order[i] << std::endl;
	} ss_graph << std::endl;

	build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 1, config.min_clip_len, ss_graph);

	std::vector<int> rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

	if (rev_topological_order.size() < n) {
		build_graph(read_seqs, order, out_edges, l_adj, l_adj_rev, 0.0, config.min_clip_len, ss_graph);

		int min_overlap = config.min_clip_len;
		for (; min_overlap <= config.read_len/2; min_overlap += 10) {
			ss_graph << "HAS A CYCLE, TRYING NO MISMATCHES AND min_overlap = " << min_overlap << std::endl;

			for (int i = 0; i < n; i++) {
				l_adj[i].erase(std::remove_if(l_adj[i].begin(), l_adj[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj[i].end());
				l_adj_rev[i].erase(std::remove_if(l_adj_rev[i].begin(), l_adj_rev[i].end(),
						[&min_overlap](edge_t& e) { return e.overlap < min_overlap; }), l_adj_rev[i].end());
				out_edges[i] = l_adj[i].size();
			}

			rev_topological_order = find_rev_topological_order(n, out_edges, l_adj_rev);

			if (rev_topological_order.size() == n) {
				ss_graph << "CYCLE DISAPPEARED!" << std::endl;
				break;
			}
		}
		if (rev_topological_order.size() < n) {
			return {"HAS_CYCLE"};
		}
	}

	std::vector<std::string> assembled_sequences;
	std::vector<bool> used(n);
	while (true) {
		// compute longest paths
		std::vector<int> best_scores(n);
		std::vector<edge_t> best_edges(n);
		for (int i : rev_topological_order) {
			if (used[i]) continue;
			if (best_scores[i] == 0 && !path_permissions[i].can_end_path) continue; // sink and cannot end path => discard
			for (edge_t& e : l_adj_rev[i]) {
				if (best_scores[e.next] < e.score + best_scores[i]) {
					best_scores[e.next] = e.score + best_scores[i];
					best_edges[e.next] = {i, e.score, e.overlap};
				}
			}
		}

		int best_score = 0, curr_vertex = 0;
		for (int i = 0; i < best_scores.size(); i++) {
			if (path_permissions[i].can_start_path && best_score < best_scores[i]) {
				best_score = best_scores[i];
				curr_vertex = i;
			}
		}
		if (best_score == 0) break;

		std::string assembled_sequence = read_seqs[curr_vertex];
		std::vector<std::string> used_reads; // track reads used to build this contig, so that we can use them for correction
		used_reads.push_back(read_seqs[curr_vertex]);
		while (best_edges[curr_vertex].overlap) {
			ss_graph << curr_vertex << " -> " << best_edges[curr_vertex].next << " , " << best_edges[curr_vertex].overlap << " ";
			ss_graph << read_seqs[curr_vertex] << std::endl;
			used[curr_vertex] = true;
			int overlap = best_edges[curr_vertex].overlap;
			curr_vertex = best_edges[curr_vertex].next;
			assembled_sequence += read_seqs[curr_vertex].substr(overlap);
			used_reads.push_back(read_seqs[curr_vertex]);
		}
		used[curr_vertex] = true;
		correct_contig(assembled_sequence, used_reads, harsh_aligner);
		assembled_sequences.push_back(assembled_sequence);
		ss_graph << "-> " << read_seqs[curr_vertex] << std::endl;
		ss_graph << assembled_sequence << std::endl << std::endl;
	}

	return assembled_sequences;
}

std::pair<StripedSmithWaterman::Alignment, StripedSmithWaterman::Alignment> remap_assembled_sequence(
		std::vector<std::string>& assembled_sequences, char* ref, int ref_len, StripedSmithWaterman::Aligner& aligner_to_base,
		bool& good_left_anchor, bool& good_right_anchor, std::string& ins_seq, std::string& full_assembled_seq) {
	std::string primary_assembled_sequence = assembled_sequences[0];

	StripedSmithWaterman::Filter filter;
	auto is_fully_aln = [](const StripedSmithWaterman::Alignment& aln, size_t seq_len) {
		if (aln.query_begin >= config.min_clip_len || aln.query_end <= seq_len-1-config.min_clip_len) return false;
		int tot_ins_sum = 0;
		for (uint32_t c : aln.cigar) {
			if (cigar_int_to_op(c) == 'I') {
				if (cigar_int_to_len(c) >= config.min_insertion_size) return false;
				tot_ins_sum += cigar_int_to_len(c);
			}
		}
		if (tot_ins_sum >= config.min_insertion_size) return false;
		return true;
	};

	int a = 0, b = primary_assembled_sequence.length();
	StripedSmithWaterman::Alignment lh_aln, rh_aln;
	std::string	assembled_seq_lh, assembled_seq_rh;
	std::string padded_assembled_seq_lh, padded_assembled_seq_rh;
	while (true) {
		int divide = (a+b)/2;
		if (divide == a) break;

		assembled_seq_lh = primary_assembled_sequence.substr(0, divide);
		assembled_seq_rh = primary_assembled_sequence.substr(divide);
		padded_assembled_seq_lh = std::string(config.clip_penalty, 'N') + assembled_seq_lh;
		padded_assembled_seq_rh = assembled_seq_rh + std::string(config.clip_penalty, 'N');
		aligner_to_base.Align(padded_assembled_seq_lh.c_str(), ref, ref_len, filter, &lh_aln, 0);
		aligner_to_base.Align(padded_assembled_seq_rh.c_str(), ref, ref_len, filter, &rh_aln, 0);

		bool lh_fully_aln = is_fully_aln(lh_aln, padded_assembled_seq_lh.length());
		bool rh_fully_aln = is_fully_aln(rh_aln, padded_assembled_seq_rh.length());

		if (lh_fully_aln && rh_fully_aln) break;
		else if (lh_fully_aln) {
			a = divide;
		} else if (rh_fully_aln) {
			b = divide;
		} else break;
	}

	good_left_anchor = lh_aln.query_begin == 0 && lh_aln.query_end-lh_aln.query_begin >= config.min_clip_len+config.clip_penalty;
	good_right_anchor = rh_aln.query_end == padded_assembled_seq_rh.length()-1 &&
			rh_aln.query_end-rh_aln.query_begin >= config.min_clip_len+config.clip_penalty;


	int lh_ref_begin = lh_aln.ref_begin + config.clip_penalty, rh_ref_end = rh_aln.ref_end - config.clip_penalty;
	int overlap_begin = std::max(lh_ref_begin, rh_aln.ref_begin),
		overlap_end = std::min(lh_aln.ref_end, rh_ref_end);
	if (good_left_anchor && good_right_anchor) {
		if (overlap_end >= overlap_begin && rh_ref_end-lh_ref_begin >= 2*config.min_clip_len) { // alignments overlap, resolve that
			std::vector<uint32_t> lh_cigar = lh_aln.cigar;
			// we want to ignore the first config.clip_penalty Xs, since they are due to the padding we added
			// however, in some cases, there are more than config.clip_penalty Xs in the first cigar op, so we need to split it
			if (cigar_int_to_op(lh_cigar[0]) == 'X' && cigar_int_to_len(lh_cigar[0]) > config.clip_penalty) {
				lh_cigar[0] = to_cigar_int(cigar_int_to_len(lh_cigar[0])-config.clip_penalty, 'X');
				lh_cigar.insert(lh_cigar.begin(), to_cigar_int(config.clip_penalty, 'X'));
			}
			std::vector<int> lh_prefix_scores = ssw_cigar_to_prefix_ref_scores(lh_cigar.data()+1, lh_cigar.size()-1);
			std::vector<uint32_t> rev_rh_cigar(rh_aln.cigar.rbegin(), rh_aln.cigar.rend());
			if (cigar_int_to_op(rev_rh_cigar[0]) == 'X' && cigar_int_to_len(rev_rh_cigar[0]) > config.clip_penalty) {
				rev_rh_cigar[0] = to_cigar_int(cigar_int_to_len(rev_rh_cigar[0])-config.clip_penalty, 'X');
				rev_rh_cigar.insert(rev_rh_cigar.begin(), to_cigar_int(config.clip_penalty, 'X'));
			}
			std::vector<int> rh_suffix_scores = ssw_cigar_to_prefix_ref_scores(rev_rh_cigar.data()+1, rev_rh_cigar.size()-1);

			int best_i = overlap_begin, best_score = INT32_MIN;
			for (int i = overlap_begin; i <= overlap_end; i++) {
				int prefix_len = i - lh_ref_begin + 1, suffix_len = rh_ref_end - i;
				if (lh_prefix_scores[prefix_len] + rh_suffix_scores[suffix_len] > best_score) {
					best_score = lh_prefix_scores[prefix_len] + rh_suffix_scores[suffix_len];
					best_i = i;
				}
			}

			aligner_to_base.Align(padded_assembled_seq_lh.c_str(), ref, best_i+1, filter, &lh_aln, 0);
			aligner_to_base.Align(padded_assembled_seq_rh.c_str(), ref+best_i+1, ref_len-(best_i+1), filter, &rh_aln, 0);
			rh_aln.ref_begin += best_i+1;
			rh_aln.ref_end += best_i+1;
		} else if (rh_aln.ref_begin <= lh_aln.ref_end) { // outwardly mapping, force one of the alignments to align to a consistent location
			if (lh_aln.sw_score >= rh_aln.sw_score) {
				aligner_to_base.Align(padded_assembled_seq_rh.c_str(), ref+lh_aln.ref_end, ref_len-lh_aln.ref_end, filter, &rh_aln, 0);
				rh_aln.ref_begin += lh_aln.ref_end+1;
				rh_aln.ref_end += lh_aln.ref_end+1;
			} else if (lh_aln.sw_score < rh_aln.sw_score) {
				aligner_to_base.Align(padded_assembled_seq_lh.c_str(), ref, rh_aln.ref_begin, filter, &lh_aln, 0);
			}
		}

		good_left_anchor = lh_aln.query_begin == 0 && lh_aln.query_end-lh_aln.query_begin >= config.min_clip_len+config.clip_penalty;
		good_right_anchor = rh_aln.query_end == padded_assembled_seq_rh.length()-1 &&
				rh_aln.query_end-rh_aln.query_begin >= config.min_clip_len+config.clip_penalty;
	}

	// if one of the anchors is bad but the other is good, we may have two half assemblies
	if (!good_left_anchor && good_right_anchor && assembled_sequences.size() >= 2) {
		padded_assembled_seq_lh = std::string(config.clip_penalty, 'N') + assembled_sequences[1];
		aligner_to_base.Align(padded_assembled_seq_lh.c_str(), ref, rh_aln.ref_begin, filter, &lh_aln, 0);
		std::string left_seq = padded_assembled_seq_lh.substr(lh_aln.query_end+1);
		if (left_seq.length() >= config.min_clip_len) {
			good_left_anchor = lh_aln.query_begin == 0 && lh_aln.query_end-lh_aln.query_begin >= config.min_clip_len+config.clip_penalty;
		}
		ins_seq = left_seq + "-" + assembled_seq_lh + padded_assembled_seq_rh.substr(0, rh_aln.query_begin);
		full_assembled_seq = assembled_sequences[1] + "-" + primary_assembled_sequence;
	} else if (good_left_anchor && !good_right_anchor && assembled_sequences.size() >= 2) {
		padded_assembled_seq_rh = assembled_sequences[1] + std::string(config.clip_penalty, 'N');
		aligner_to_base.Align(padded_assembled_seq_rh.c_str(), ref+lh_aln.ref_end, ref_len-lh_aln.ref_end, filter, &rh_aln, 0);
		rh_aln.ref_begin += lh_aln.ref_end;
		rh_aln.ref_end += lh_aln.ref_end;
		std::string right_seq = padded_assembled_seq_rh.substr(0, rh_aln.query_begin);
		if (right_seq.length() >= config.min_clip_len) {
			good_right_anchor = rh_aln.query_end == padded_assembled_seq_rh.length()-1 &&
					rh_aln.query_end-rh_aln.query_begin >= config.min_clip_len+config.clip_penalty;
		}
		ins_seq = padded_assembled_seq_lh.substr(lh_aln.query_end+1) + assembled_seq_rh + "-" + right_seq;
		full_assembled_seq = primary_assembled_sequence + "-" + assembled_sequences[1];
	} else {
		ins_seq = padded_assembled_seq_lh.substr(lh_aln.query_end+1) + padded_assembled_seq_rh.substr(0, rh_aln.query_begin);
		full_assembled_seq = primary_assembled_sequence;
	}

	return {lh_aln, rh_aln};
}

std::vector<std::string> assemble_sequences(std::string contig_name, reads_cluster_t* r_cluster, reads_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, StripedSmithWaterman::Aligner& harsh_aligner) {
	std::vector<seq_w_pp_t> left_stable_read_seqs, unstable_read_seqs, right_stable_read_seqs;
	std::unordered_set<std::string> used_ls, used_us, used_rs;
	if (r_cluster->clip_cluster) left_stable_read_seqs.push_back({r_cluster->clip_cluster->full_seq, true, false});
	for (bam1_t* read : r_cluster->reads) {
		std::string seq = get_sequence(read);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		rc(mate_seq);
		if (!used_ls.count(seq)) {
			left_stable_read_seqs.push_back({seq, !is_left_clipped(read, 0), true});
			used_ls.insert(seq);
		}
		if (!used_us.count(mate_seq)) {
			unstable_read_seqs.push_back({mate_seq, true, true});
			used_us.insert(mate_seq);
		}
	}
	for (bam1_t* read : l_cluster->reads) {
		std::string seq = get_sequence(read);
		std::string mate_seq = get_mate_seq(read, mateseqs);
		if (!used_rs.count(seq)) {
			right_stable_read_seqs.push_back({seq, true, !is_right_clipped(read, 0)});
			used_rs.insert(seq);
		}
		if (!used_us.count(mate_seq)) {
			unstable_read_seqs.push_back({mate_seq, true, true});
			used_us.insert(mate_seq);
		}
	}
	if (l_cluster->clip_cluster) right_stable_read_seqs.push_back({l_cluster->clip_cluster->full_seq, false, true});

	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		if (!used_us.count(read_seq)) {
			unstable_read_seqs.push_back({read_seq, !is_left_clipped(read, 0), !is_right_clipped(read, 0)});
			used_us.insert(read_seq);
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		std::string read_seq = get_sequence(read, true);
		if (!used_us.count(read_seq)) {
			unstable_read_seqs.push_back({read_seq, !is_left_clipped(read, 0), !is_right_clipped(read, 0)});
			used_us.insert(read_seq);
		}
	}

	if (unstable_read_seqs.size() + left_stable_read_seqs.size() + right_stable_read_seqs.size() >= TOO_MANY_READS) return {"TOO_MANY_READS"};

	std::stringstream ss_graph;
	std::vector<std::string> assembled_sequences = assemble_reads(left_stable_read_seqs, unstable_read_seqs, right_stable_read_seqs,
			harsh_aligner, config, ss_graph);

	return assembled_sequences;
}

insertion_t* assemble_insertion(std::string& contig_name, chr_seqs_map_t& contigs,
		reads_cluster_t* r_cluster, reads_cluster_t* l_cluster,
		std::unordered_map<std::string, std::string>& mateseqs, std::unordered_map<std::string, std::string>& matequals,
		StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
		std::vector<bam1_t*>& assembled_reads) {

	std::vector<std::string> assembled_sequences = assemble_sequences(contig_name, r_cluster, l_cluster, mateseqs, harsh_aligner);

	std::string ins_full_id = "NO_ID";
	if (assembled_sequences.empty()) {
		failed_assembly_mtx.lock();
		assembly_failed_no_seq << ins_full_id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_no_seq << contig_name << " " << l_cluster->start() << " - INS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "HAS_CYCLE") {
		failed_assembly_mtx.lock();
		assembly_failed_cycle_writer << ins_full_id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_cycle_writer << contig_name << " " << l_cluster->start() << " - INS HAS_CYCLE" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	} else if (assembled_sequences[0] == "TOO_MANY_READS") {
		failed_assembly_mtx.lock();
		assembly_failed_too_many_reads_writer << ins_full_id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_too_many_reads_writer << contig_name << " " << l_cluster->start() << " - INS TOO_MANY_READS" << std::endl;
		failed_assembly_mtx.unlock();
		return NULL;
	}

	// remap assembled sequences to find breakpoints
	char* contig = contigs.get_seq(contig_name);
	int contig_len = contigs.get_len(contig_name);

	std::string assembled_sequence = assembled_sequences[0];

	// divide the assembled sequence into two and remap them
	int extend = 2*(assembled_sequence.length() + config.clip_penalty);
	int remap_region_start = std::min(r_cluster->end(), l_cluster->start())-extend;
	remap_region_start = std::max(0, remap_region_start);
	int remap_region_end = std::max(r_cluster->end(), l_cluster->start())+extend;
	remap_region_end = std::min(remap_region_end, contig_len-1);
	int remap_region_len = remap_region_end - remap_region_start;

	std::string ins_seq, full_assembled_seq;
	bool good_left_anchor, good_right_anchor;
	auto alns = remap_assembled_sequence(assembled_sequences, contig+remap_region_start, remap_region_end-remap_region_start, aligner_to_base,
			good_left_anchor, good_right_anchor, ins_seq, full_assembled_seq);
	StripedSmithWaterman::Alignment lh_aln = alns.first, rh_aln = alns.second;
	if (ins_seq.find("N") != std::string::npos) return NULL;

	std::string mh_seq;
	int ins_start = remap_region_start + lh_aln.ref_end, ins_end = remap_region_start + rh_aln.ref_begin - 1;
	if (ins_start > ins_end) {
		int mh_len = ins_start - ins_end;
		char* mh_seq_cstr = new char[mh_len+1];
		strncpy(mh_seq_cstr, contigs.get_seq(contig_name)+ins_end, mh_len);
		mh_seq_cstr[mh_len] = '\0';
		for (int i = 0; i < mh_len; i++) {
			mh_seq_cstr[i] = toupper(mh_seq_cstr[i]);
		}
		mh_seq = mh_seq_cstr;
		ins_start = ins_end;
		delete[] mh_seq_cstr;
	}
	std::string ins_seq_w_mh = mh_seq + ins_seq;

	insertion_t* ins = new insertion_t(contig_name, ins_start, ins_end, 0, 0, 0, 0, 0, 0, 0, ins_seq_w_mh);
//	ins->id = ins_full_id;

	if (!good_left_anchor || !good_right_anchor) {
		assembly_failed_bad_anchors_writer << ins->id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_bad_anchors_writer << contig_name << " " << l_cluster->start() << " - INS " << assembled_sequence << " ";
		assembly_failed_bad_anchors_writer << contig_name << ":" << remap_region_start << "-" << remap_region_end << " ";
		assembly_failed_bad_anchors_writer << lh_aln.cigar_string << " " << rh_aln.cigar_string << " ";
		assembly_failed_bad_anchors_writer << std::endl;
		return NULL;
	} else if (mh_seq.length() > config.max_mh_len) {
		assembly_failed_mh_too_long << ins->id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_mh_too_long << contig_name << " " << l_cluster->start() << " - INS " << full_assembled_seq << " ";
		assembly_failed_mh_too_long << remap_region_start + lh_aln.ref_begin << " " << remap_region_start + lh_aln.ref_end << " " << lh_aln.cigar_string << " ";
		assembly_failed_mh_too_long << remap_region_start + rh_aln.ref_begin << " " << remap_region_start + rh_aln.ref_end << " " << rh_aln.cigar_string << " ";
		assembly_failed_mh_too_long << contig_name << ":" << remap_region_start << "-" << remap_region_end << std::endl;
		return NULL;
	} else if (ins->ins_seq.length() < config.min_insertion_size) {
		assembly_failed_lt50bp << ins->id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_failed_lt50bp << contig_name << " " << l_cluster->start() << " - INS " << ins->ins_seq << std::endl;
		return NULL;
	} else {
		assembly_succeeded << ins->id << " " << contig_name << " " << r_cluster->end() << " + ";
		assembly_succeeded << contig_name << " " << l_cluster->start() << " - INS " << assembled_sequence << std::endl;
	}

	StripedSmithWaterman::Filter filter;
	if (r_cluster->clip_cluster) {
		StripedSmithWaterman::Alignment aln;
		std::string padded_clip_fullseq = std::string(config.clip_penalty, 'N') + r_cluster->clip_cluster->full_seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_clip_fullseq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.max_seq_error)) {
			ins->rc_fwd_reads = r_cluster->clip_cluster->c->a1.fwd_sc_reads;
			ins->rc_rev_reads = r_cluster->clip_cluster->c->a1.rev_sc_reads;
		}
	}
	for (bam1_t* read : r_cluster->reads) {
		std::string mate_seq = get_mate_seq(read, mateseqs);
		std::string mate_qual = get_mate_qual(read, matequals);
		rc(mate_seq);
		mate_qual = std::string(mate_qual.rbegin(), mate_qual.rend());
		mate_seq = std::string(config.clip_penalty, 'N') + mate_seq + std::string(config.clip_penalty, 'N');
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(mate_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.max_seq_error, mate_qual, stats.get_min_avg_base_qual())) {
			ins->r_disc_pairs++;
			ins->rc_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	ins->rc_avg_nm /= ins->r_disc_pairs;
	for (bam1_t* read : l_cluster->reads) {
		std::string mate_seq = std::string(config.clip_penalty, 'N') + get_mate_seq(read, mateseqs) + std::string(config.clip_penalty, 'N');
		std::string mate_qual = get_mate_qual(read, matequals);
		StripedSmithWaterman::Alignment aln;
		harsh_aligner.Align(mate_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.max_seq_error, mate_qual, stats.get_min_avg_base_qual())) {
			ins->l_disc_pairs++;
			ins->lc_avg_nm += bam_aux2i(bam_aux_get(read, "NM"));
			bam1_t* d = bam_dup1(read);
			d->core.mpos = 0;
			d->core.mtid = d->core.tid;
			d->core.flag |= BAM_FMUNMAP;
			assembled_reads.push_back(d);
		}
	}
	ins->lc_avg_nm /= ins->l_disc_pairs;
	if (l_cluster->clip_cluster) {
		StripedSmithWaterman::Alignment aln;
		std::string padded_clip_fullseq = std::string(config.clip_penalty, 'N') + l_cluster->clip_cluster->full_seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_clip_fullseq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		if (accept(aln, config.max_seq_error)) {
			ins->lc_fwd_reads = l_cluster->clip_cluster->c->a1.fwd_sc_reads;
			ins->lc_rev_reads = l_cluster->clip_cluster->c->a1.rev_sc_reads;
		}
	}

	// start and end of inserted sequence within the full assembled sequence
	int ins_seq_start = lh_aln.query_end+1-config.clip_penalty, ins_seq_end = ins_seq_start + ins_seq.length();
	for (bam1_t* read : r_cluster->semi_mapped_reads) {
		StripedSmithWaterman::Alignment aln;
		std::string read_seq = get_sequence(read, true);
		rc(read_seq);
		std::string padded_seq = std::string(config.clip_penalty, 'N') + read_seq + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		std::string qual_ascii = get_qual_ascii(read, true);
		qual_ascii = std::string(qual_ascii.rbegin(), qual_ascii.rend());
		int aln_ref_begin = aln.ref_begin + (is_left_clipped(aln) ? 0 : config.clip_penalty);
		int aln_ref_end = aln.ref_end - (is_right_clipped(aln) ? 0 : config.clip_penalty);
		if (overlap(ins_seq_start, ins_seq_end, aln_ref_begin, aln_ref_end) >= config.min_clip_len
				&& accept(aln, config.max_seq_error, qual_ascii, stats.get_min_avg_base_qual())) {
			ins->r_disc_pairs++;
			assembled_reads.push_back(bam_dup1(read));
		}
	}
	for (bam1_t* read : l_cluster->semi_mapped_reads) {
		StripedSmithWaterman::Alignment aln;
		std::string padded_seq = std::string(config.clip_penalty, 'N') + get_sequence(read, true) + std::string(config.clip_penalty, 'N');
		harsh_aligner.Align(padded_seq.c_str(), full_assembled_seq.c_str(), full_assembled_seq.length(), filter, &aln, 0);
		std::string qual_ascii = get_qual_ascii(read, true);
		int aln_ref_begin = aln.ref_begin + (is_left_clipped(aln) ? 0 : config.clip_penalty);
		int aln_ref_end = aln.ref_end - (is_right_clipped(aln) ? 0 : config.clip_penalty);
		if (overlap(ins_seq_start, ins_seq_end, aln_ref_begin, aln_ref_end) >= config.min_clip_len
				&& accept(aln, config.max_seq_error, qual_ascii, stats.get_min_avg_base_qual())) {
			ins->l_disc_pairs++;
			assembled_reads.push_back(bam_dup1(read));
		}
	}

	for (bam1_t* read : assembled_reads) {
		bam_aux_update_str(read, "ID", ins->id.length(), ins->id.c_str());
	}



	ins->left_anchor = std::to_string(remap_region_start + lh_aln.ref_begin) + "-" + std::to_string(remap_region_start + lh_aln.ref_end);
	ins->right_anchor = std::to_string(remap_region_start + rh_aln.ref_begin) + "-" + std::to_string(remap_region_start + rh_aln.ref_end);
	return ins;
}


#endif /* ASSEMBLE_H_ */
