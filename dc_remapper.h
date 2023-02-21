#ifndef DC_REMAPPER_H_
#define DC_REMAPPER_H_

#include "cluster.h"


struct clip_cluster_t {
    cluster_t* c;
    std::string clipped_seq, full_seq;

    clip_cluster_t(cluster_t* c, std::string& clipped_seq, std::string& full_seq) : c(c), clipped_seq(clipped_seq), full_seq(full_seq) {}

    int pos() {
    	return c->a1.pos();
    }

    std::string to_string() {
        return std::to_string(c->a1.start) + "-" + std::to_string(c->a1.end) + "," + std::to_string(c->a1.sc_reads());
    }
};

struct reads_cluster_t {
    std::vector<bam1_t*> reads;
    std::vector<bam1_t*> semi_mapped_reads;
    clip_cluster_t* clip_cluster = NULL;
    bool used = false;

    void add_read(bam1_t* read) {
        reads.push_back(read);
    }
    void add_semi_mapped_reads(bam1_t* read) {
    	semi_mapped_reads.push_back(read);
    }
    void add_clip_cluster(clip_cluster_t* clip_cluster) {
        if (!clip_cluster) return;
        if (!this->clip_cluster || this->clip_cluster->c->a1.sc_reads() < clip_cluster->c->a1.sc_reads()) {
            this->clip_cluster = clip_cluster;
        }
    }

    hts_pos_t start() {
        if (clip_cluster && clip_cluster->c->a1.dir == 'L') {
            return clip_cluster->c->a1.start;
        } else {
            hts_pos_t _start = INT32_MAX;
            for (bam1_t* read : reads) {
                _start = std::min(_start, read->core.pos);
            }
            for (bam1_t* read : semi_mapped_reads) {
				_start = std::min(_start, read->core.pos);
			}
            if (clip_cluster && clip_cluster->c->a1.start < _start) _start = clip_cluster->c->a1.start;
            return _start;
        }
    }
    hts_pos_t end() {
        if (clip_cluster && clip_cluster->c->a1.dir == 'R') {
            return clip_cluster->c->a1.end;
        } else {
            hts_pos_t _end = 0;
            for (bam1_t* read : reads) {
                _end = std::max(_end, bam_endpos(read));
            }
            for (bam1_t* read : semi_mapped_reads) {
				_end = std::max(_end, bam_endpos(read));
			}
            if (clip_cluster && clip_cluster->c->a1.end > _end) _end = clip_cluster->c->a1.end;
            return _end;
        }
    }
    char dir() {
        if (reads.empty()) return 'N';
        else return bam_is_rev(reads[0]) ? 'L' : 'R';
    }

    std::string to_string() {
        std::string str = std::to_string(start()) + " " + std::to_string(end()) + " ";
        if (clip_cluster) str += clip_cluster->to_string();
        return str;
    }

    void deallocate_reads() {
        for (bam1_t* read : reads) {
            bam_destroy1(read);
        }
        for (bam1_t* read : semi_mapped_reads) {
        	bam_destroy1(read);
        }
        reads.clear();
        semi_mapped_reads.clear();
    }

    bool empty() {
    	return reads.empty() && !clip_cluster;
    }
};

std::string get_mate_seq(bam1_t* read,  std::unordered_map<std::string, std::string>& mateseqs) {
    std::string qname = bam_get_qname(read);
    if (is_samechr(read)) {
        if (read->core.flag & BAM_FREAD1) qname += "_2";
        else qname += "_1";
    }
    if (!mateseqs.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return mateseqs[qname];
}
std::string get_mate_qual(bam1_t* read,  std::unordered_map<std::string, std::string>& matequals) {
    std::string qname = bam_get_qname(read);
    if (is_samechr(read)) {
        if (read->core.flag & BAM_FREAD1) qname += "_2";
        else qname += "_1";
    }
    if (!matequals.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return matequals[qname];
}


#endif /* DC_REMAPPER_H_ */
