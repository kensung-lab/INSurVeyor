#ifndef SMALLINSFINDER_UTILS_H
#define SMALLINSFINDER_UTILS_H

#include <unordered_map>
#include <vector>
#include <fstream>
#include <unistd.h>

#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read)

struct config_t {
    int threads;
    int min_insertion_size, max_insertion_size;
    int max_clipped_pos_dist, min_clip_len, min_stable_mapq;
    double max_seq_error;

    int clip_penalty = 7;
    int max_mh_len = 100;

    int max_is, read_len;

    void parse(std::string config_fname) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_fname);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = stoi(config_params["threads"]);
        min_insertion_size = stoi(config_params["min_insertion_size"]);
        max_insertion_size = stoi(config_params["max_insertion_size"]);
        max_clipped_pos_dist = stoi(config_params["max_clipped_pos_dist"]);
        min_clip_len = stoi(config_params["min_clip_len"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        min_stable_mapq = stoi(config_params["min_stable_mapq"]);

        max_is = stoi(config_params["max_is"]);
        read_len = stoi(config_params["read_len"]);
    };
};

struct stats_t {
    int min_depth, median_depth, max_depth, min_avg_base_qual;

    void parse_stats(std::string stats_fname) {
		std::unordered_map<std::string, std::string> params;
		std::ifstream fin(stats_fname);
		std::string name, value;
		while (fin >> name >> value) {
			params[name] = value;
		}
		fin.close();

		min_depth = stoi(params["min_depth"]);
		median_depth = stoi(params["median_depth"]);
		max_depth = stoi(params["max_depth"]);
		min_avg_base_qual = stoi(params["min_avg_base_qual"]);
	}
};


struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    void parse(std::string workdir) {
    	std::ifstream fin(workdir + "/contig_map");
		std::string name;
		int id = 0;
		while (fin >> name) {
			name_to_id[name] = id;
			id_to_name.push_back(name);
			id++;
		}
    }

    size_t size() {return id_to_name.size();}
    std::string get_name(size_t id) {return id_to_name[id];};
    size_t get_id(std::string& name) {return name_to_id[name];};
};

struct chr_seq_t {
    char* seq;
    hts_pos_t len;

    chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};

struct suffix_prefix_aln_t {
    int overlap, score, mismatches;

    suffix_prefix_aln_t(int overlap, int score, int mismatches) : overlap(overlap), score(score), mismatches(mismatches) {}
};

// Finds the best alignment between a suffix of s1 and a prefix of s2
// Disallows gaps
suffix_prefix_aln_t aln_suffix_prefix(std::string& s1, std::string& s2, int match_score, int mismatch_score, double max_seq_error,
                                      int min_overlap = 1, int max_overlap = INT32_MAX, int max_mismatches = INT32_MAX) {
    int best_score = 0, best_aln_mismatches = 0;
    int overlap = 0;

    for (int i = std::max(0, (int) s1.length()-max_overlap); i < s1.length()-min_overlap+1; i++) {
        if (i+s2.length() < s1.length()) continue;

        int sp_len = s1.length()-i;
        if (best_score >= sp_len*match_score) break; // current best score is unbeatable

        int mismatches = 0;
        const char* s1_suffix = s1.data()+i;
        const char* s2_prefix = s2.data();
        while (*s1_suffix) {
            if (*s1_suffix != *s2_prefix) mismatches++;
            s1_suffix++; s2_prefix++;
        }

        int score = (sp_len-mismatches)*match_score + mismatches*mismatch_score;

        int max_acceptable_mm = max_seq_error == 0.0 ? 0 : std::max(1.0, sp_len*max_seq_error);
        if (best_score < score && mismatches <= max_acceptable_mm && mismatches <= max_mismatches) {
            best_score = score;
            best_aln_mismatches = mismatches;
            overlap = sp_len;
        }
    }
    return suffix_prefix_aln_t(overlap, best_score, best_aln_mismatches);
}

template<typename T>
inline T max(T a, T b, T c) { return std::max(std::max(a,b), c); }

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

int max_pos(int* a, int n) {
	int pos = 0;
	for (int i = 1; i < n; i++) {
		if (a[i] > a[pos]) pos = i;
	}
	return pos;
}

bool file_exists(std::string& fname) {
	return std::ifstream(fname).good();
}

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

struct insertion_t {
    std::string id;
    std::string chr;
    hts_pos_t start, end;
    int overlap, r_disc_pairs, l_disc_pairs, rc_reads, lc_reads;
    int r_conc_pairs = 0, l_conc_pairs = 0, median_lf_cov = 0, median_rf_cov = 0;
    int mh_len = 0;
    std::string ins_seq;
    std::string left_anchor, right_anchor, left_anchor_cigar, right_anchor_cigar;
    bool left_bp_precise = false, right_bp_precise = false;
    double rc_avg_nm = 0.0, lc_avg_nm = 0.0;

    insertion_t(std::string chr, hts_pos_t start, hts_pos_t end, int r_disc_pairs, int l_disc_pairs,
    		int rc_reads, int lc_reads, int overlap, std::string ins_seq) :
    chr(chr), start(start), end(end), r_disc_pairs(r_disc_pairs), l_disc_pairs(l_disc_pairs), rc_reads(rc_reads), lc_reads(lc_reads),
	overlap(overlap), ins_seq(ins_seq) {}
};


int score(char a, char b, int match_score, int mismatch_penalty) {
	return (toupper(a) == toupper(b) || a == 'N' || b == 'N') ? match_score : mismatch_penalty;
}
int* smith_waterman_gotoh(const char* ref, int ref_len, const char* read, int read_len, int match_score, int mismatch_penalty, int gap_open, int gap_extend) {
	const int INF = 1000000;

	int** dab = new int*[ref_len+1];
	int** dag = new int*[ref_len+1];
	int** dgb = new int*[ref_len+1];
	for (int i = 0; i <= ref_len; i++) {
		dab[i] = new int[read_len+1];
		dag[i] = new int[read_len+1];
		dgb[i] = new int[read_len+1];
		std::fill(dab[i], dab[i]+read_len+1, 0);
		std::fill(dag[i], dag[i]+read_len+1, 0);
		std::fill(dgb[i], dgb[i]+read_len+1, 0);
	}

	for (int i = 1; i <= ref_len; i++) {
		dab[i][0] = -INF;
		dag[i][0] = -INF;
		dgb[i][0] = gap_open + (i-1)*gap_extend;
	}
	for (int i = 1; i <= read_len; i++) {
		dab[0][i] = -INF;
		dag[0][i] = gap_open + (i-1)*gap_extend;
		dgb[0][i] = -INF;
	}

	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			dab[i][j] = score(ref[i-1], read[j-1], match_score, mismatch_penalty) + max(dab[i-1][j-1], dag[i-1][j-1], dgb[i-1][j-1], 0);
			dag[i][j] = max(gap_open + dab[i][j-1], gap_extend + dag[i][j-1], gap_open + dgb[i][j-1]);
			dgb[i][j] = max(gap_open + dab[i-1][j], gap_open + dag[i-1][j], gap_extend + dgb[i-1][j]);
		}
	}

	int* prefix_scores = new int[read_len];
	std::fill(prefix_scores, prefix_scores+read_len, 0);
	for (int i = 1; i <= ref_len; i++) {
		for (int j = 1; j <= read_len; j++) {
			prefix_scores[j-1] = std::max(prefix_scores[j-1], dab[i][j]);
		}
	}

	for (int i = 0; i <= ref_len; i++) {
		delete[] dab[i];
		delete[] dag[i];
		delete[] dgb[i];
	}
	delete[] dab;
	delete[] dag;
	delete[] dgb;

	for (int i = 1; i < read_len; i++) {
		prefix_scores[i] = std::max(prefix_scores[i], prefix_scores[i-1]);
	}

	return prefix_scores;
}

int get_left_clip_size(const StripedSmithWaterman::Alignment& aln) {
    uint32_t l = aln.cigar[0];
    return cigar_int_to_op(l) == 'S' ? cigar_int_to_len(l) : 0;
}
int get_right_clip_size(const StripedSmithWaterman::Alignment& aln) {
    uint32_t r = aln.cigar[aln.cigar.size()-1];
    return cigar_int_to_op(r) == 'S' ? cigar_int_to_len(r) : 0;
}


#endif //SMALLINSFINDER_UTILS_H
