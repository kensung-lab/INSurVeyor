#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include <iostream>
#include <atomic>

#include "sam_utils.h"
#include "utils.h"
#include "htslib/sam.h"

struct anchor_t {
    static constexpr const char* pattern = "%c:%d:%d:%d:%d";

    char dir;
    int contig_id;
    int start, end;
    int sc_reads;

    anchor_t() {}
    anchor_t(char dir, int contig_id, int start, int end, int sc_reads) : dir(dir), contig_id(contig_id), start(start),
                                                                          end(end), sc_reads(sc_reads) {}
    anchor_t(char* s) {
        sscanf(s, pattern, &dir, &contig_id, &start, &end, &sc_reads);
    }

    int pos() { return dir == 'L' ? start : end; }

    static bool check_clipped(anchor_t& clipped, anchor_t& other, config_t& config) {
        if (clipped.dir == 'L') {
            return clipped.pos() <= other.pos()+config.max_clipped_pos_dist;
        } else {
            return clipped.pos() >= other.pos()-config.max_clipped_pos_dist;
        }
    }

    static bool can_merge(anchor_t& a1, anchor_t& a2, config_t& config) {
        if (a1.sc_reads > 0 && !check_clipped(a1, a2, config)) return false;
        if (a2.sc_reads > 0 && !check_clipped(a2, a1, config)) return false;
        return distance(a1, a2) <= config.max_is;
    }

    static anchor_t merge(anchor_t& a1, anchor_t& a2) {
        return anchor_t(a1.dir, a1.contig_id, std::min(a1.start, a2.start),
                        std::max(a1.end,  a2.end), a1.sc_reads+a2.sc_reads);
    }

    static int distance(anchor_t& a1, anchor_t& a2) {
        if (a1.contig_id != a2.contig_id || a1.dir != a2.dir) return INT_MAX;
        return std::max(a1.end, a2.end) - std::min(a1.start, a2.start);
    }

    int size() {
        return end-start+1;
    }
};
bool operator < (const anchor_t& a1, const anchor_t& a2) {
    if (a1.start != a2.start) return a1.start < a2.start;
    return a1.end < a2.end;
}

struct cluster_t {
    static constexpr const char* scan_pattern = "%[^-]-%[^-]-%d";

    int id; // not necessarily set, use if needed
    anchor_t a1, a2;
    int disc_pairs;
    bool dead = false;

    cluster_t(const anchor_t &a1, const anchor_t &a2, int disc_pairs) : a1(a1), a2(a2), disc_pairs(disc_pairs) {}

    cluster_t(std::string& s) {
        char anchor1[100], anchor2[100];
        sscanf(s.c_str(), scan_pattern, anchor1, anchor2, &disc_pairs);
        a1 = anchor_t(anchor1);
        a2 = anchor_t(anchor2);
    }

    cluster_t(cluster_t* c) {
        this->id = c->id;
        this->a1 = c->a1;
        this->a2 = c->a2;
        this->disc_pairs = c->disc_pairs;
        this->dead = c->dead;
    }

    static bool can_merge(cluster_t* c1, cluster_t* c2, config_t& config) {
        return anchor_t::can_merge(c1->a1, c2->a1, config) && anchor_t::can_merge(c1->a2, c2->a2, config);
    }
    static int distance(cluster_t* c1, cluster_t* c2) {
        if (c1->a1.contig_id != c2->a1.contig_id || c1->a1.dir != c2->a1.dir ||
            c1->a2.contig_id != c2->a2.contig_id || c1->a2.dir != c2->a2.dir) {
            return INT32_MAX;
        }
        return anchor_t::distance(c1->a1, c2->a1) + anchor_t::distance(c1->a2, c2->a2);
    }

    static cluster_t* merge(cluster_t* c1, cluster_t* c2) {
        return new cluster_t(anchor_t::merge(c1->a1, c2->a1), anchor_t::merge(c1->a2, c2->a2),
                           c1->disc_pairs+c2->disc_pairs);
    }
};


struct cc_distance_t {
    int distance;
    cluster_t* c1,* c2;

    cc_distance_t(int distance, cluster_t *c1, cluster_t *c2) : distance(distance), c1(c1), c2(c2) {}
};
bool operator < (const cc_distance_t& ccd1, const cc_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance > ccd2.distance;
}


struct breakpoint_t {
    static constexpr const char* pattern = "%c:%d:%d:%d:%d:%d";

    char dir;
    int contig_id, start, end;
    int sc_reads;
    int spanning_reads = 0;

    breakpoint_t() {}

    breakpoint_t(anchor_t anchor) : dir(anchor.dir), contig_id(anchor.contig_id),
                                    start(anchor.start), end(anchor.end), sc_reads(anchor.sc_reads) {}

    breakpoint_t(char* s) {
        sscanf(s, pattern, &dir, &contig_id, &start, &end, &sc_reads, &spanning_reads);
    }

    int pos() { return dir == 'L' ? start : end; }

    std::string to_string() {
        char buffer[100];
        sprintf(buffer, pattern, dir, contig_id, start, end, sc_reads, spanning_reads);
        return std::string(buffer);
    }
};

std::atomic<int> pred_id(0);

struct prediction_t {
    static constexpr const char* scan_pattern = "%d,%[^,],%[^,],%[^,],%d,%d,%d,%lf,%lf,%d";
    static constexpr const char* print_pattern = "%d,%s,%s,%s,%d,%d,%d,%lf,%lf,%d";

    int id;
    breakpoint_t bp1, bp2;
    int disc_pairs;
    int size = INT32_MAX, conf_ival = 0;
    double pval = -1.0, shift_pval = 1.0;
    int stable_depth = 0;

    prediction_t() {}

    prediction_t(cluster_t* c) : id(pred_id++), bp1(c->a1), bp2(c->a2), disc_pairs(c->disc_pairs) {}

    prediction_t(std::string& line) {
        char breakpoint1[1000], breakpoint2[1000];
        char svt[10];
        sscanf(line.c_str(), scan_pattern, &id, breakpoint1, breakpoint2, svt, &disc_pairs, &size, &conf_ival, &pval,
               &shift_pval, &stable_depth);
        if (id >= pred_id) pred_id = id+1;
        bp1 = breakpoint_t(breakpoint1);
        bp2 = breakpoint_t(breakpoint2);
    }

    std::string to_str() {
        char buffer[10000];
        sprintf(buffer, print_pattern, id, bp1.to_string().c_str(), bp2.to_string().c_str(), "TRA", disc_pairs, size,
                conf_ival, pval, shift_pval, stable_depth);
        return std::string(buffer);
    }
};

#endif //SURVEYOR_CLUSTER_H
