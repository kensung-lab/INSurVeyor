#ifndef SMALLINSFINDER_SAM_UTILS_H
#define SMALLINSFINDER_SAM_UTILS_H

#include <vector>
#include <sstream>
#include <algorithm>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "utils.h"

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}
bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}
bool is_proper_pair(bam1_t* r) {
	return is_primary(r) && (r->core.flag & BAM_FPROPER_PAIR);
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}
bool is_samestr(bam1_t* r) {
    return (r->core.flag & BAM_FREVERSE) == (r->core.flag & BAM_FMREVERSE);
}
bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000 || is_unmapped(r) != is_mate_unmapped(r);
}

bool is_left_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' && bam_cigar_oplen(cigar[0]) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' && bam_cigar_oplen(cigar[r->core.n_cigar-1]) >= min_clip_len;
}

char* get_mc_tag(bam1_t* r) {
	const uint8_t* mc_tag = bam_aux_get(r, "MC");
	if (mc_tag == NULL) {
		std::cerr << "Warning: read " << std::string(bam_get_qname(r)) << " does not have the MC tag.";
		return (char*) "";
	}
	return bam_aux2Z(mc_tag);
}

bool is_mate_left_clipped(bam1_t* r) {
    char* mc_tag_str = get_mc_tag(r);
    int i = 0;
    while (mc_tag_str[i] >= '0' && mc_tag_str[i] <= '9') i++;
    return mc_tag_str[i] == 'S';
}
bool is_mate_right_clipped(bam1_t* r) {
	char* mc_tag_str = get_mc_tag(r);
    int i = strlen(mc_tag_str)-1;
    return mc_tag_str[i] == 'S';
}
bool is_mate_clipped(bam1_t* r) {
	return is_mate_left_clipped(r) && is_mate_right_clipped(r);
}

int get_left_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]): 0;
}
int get_right_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]): 0;
}

int get_endpoint(bam1_t* r) {
    return bam_is_rev(r) ? r->core.pos : bam_endpos(r);
}

int get_mate_endpos(const bam1_t* r) {
    uint8_t *mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}

int64_t get_mq(bam1_t* r) {
    uint8_t* mq = bam_aux_get(r, "MQ");
    if (mq == NULL) {
    	if ((r->core.flag & BAM_FUNMAP) == 0 && (r->core.flag & BAM_FMUNMAP) == 0) {
			std::cerr << "Warning: read pair " << bam_get_qname(r) << " does not have an MQ tag. Please include it." << std::endl;
			std::cerr << (r->core.flag & BAM_FMUNMAP) << std::endl;
    	}
    	return 0;
    }
    return bam_aux2i(mq);
}

double avg_qual(bam1_t* read) {
	double avg_qual = 0;
	uint8_t* qual = bam_get_qual(read);
	for (int i = 0; i < read->core.l_qseq; i++) {
		avg_qual += qual[i];
	}
	return avg_qual/read->core.l_qseq;
}
double avg_qual(std::string& qual_ascii, int offset = 33) {
	double tot = 0;
	for (char c : qual_ascii) tot += int(c)-offset;
	return tot/qual_ascii.length();
}
int median_qual(std::string qual_ascii, int offset = 33) {
	std::sort(qual_ascii.begin(), qual_ascii.end());
	return int(qual_ascii[qual_ascii.size()/2])-offset;
}

std::string get_cigar_code(bam1_t* r) {
    const uint32_t* cigar = bam_get_cigar(r);
    std::stringstream ss;
    for (int i = 0; i < r->core.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

void rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
		char c = std::toupper(read[i]);
		if (c == 'A') read[i] = 'T';
		else if (c == 'C') read[i] = 'G';
		else if (c == 'G') read[i] = 'C';
		else if (c == 'T') read[i] = 'A';
		else c = 'N';
	}
}
void rc(char* read) {
    int len = strlen(read);
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
    	char c = std::toupper(read[i]);
        if (c == 'A') read[i] = 'T';
        else if (c == 'C') read[i] = 'G';
        else if (c == 'G') read[i] = 'C';
        else if (c == 'T') read[i] = 'A';
        else c = 'N';
    }
}

char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}
std::string get_sequence(bam1_t* r, bool fastq_seq = false) {
    char seq[100000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
    return std::string(seq);
}
std::string get_qual_ascii(bam1_t* r, bool fastq_seq = false) {
	uint8_t* qual = bam_get_qual(r);
	std::string qual_ascii(r->core.l_qseq, ' ');
	for (int i = 0; i < r->core.l_qseq; i++) {
		qual_ascii[i] = char(33 + qual[i]);
	}
	if (fastq_seq && bam_is_rev(r)) qual_ascii = std::string(qual_ascii.rbegin(), qual_ascii.rend());
	return qual_ascii;
}

bool is_homopolymer(const char* seq, int len) {
	int a = 0, c = 0, g = 0, t = 0;
	for (int i = 0; i < len; i++) {
		char b = std::toupper(seq[i]);
		if (b == 'A') a++;
		else if (b == 'C') c++;
		else if (b == 'G') g++;
		else if (b == 'T') t++;
	}
	return max(a, c, g, t)/double(a+c+g+t) >= 0.8;
}

bool is_homopolymer(std::string seq) {
	return is_homopolymer(seq.data(), seq.length());
}

samFile* open_writer(std::string path, bam_hdr_t* header) {
    samFile* remapped_file = sam_open(path.c_str(), "wb");
    if (remapped_file == NULL) {
        throw "Unable to open " + path;
    }
    if (sam_hdr_write(remapped_file, header) != 0) {
        throw "Could not write file " + std::string(remapped_file->fn);
    }
    return remapped_file;
}

void write_and_index_file(std::vector<bam1_t*>& reads, std::string path, bam_hdr_t* header) {
    samFile* file = open_writer(path, header);
    if (file == NULL) {
        throw "Unable to open " + path;
    }

    // write reads
    std::sort(reads.begin(), reads.end(), [](bam1_t *r1, bam1_t *r2) { return r1->core.pos < r2->core.pos; });
    for (bam1_t* r : reads) {
        int ok = sam_write1(file, header, r);
        if (ok < 0) throw "Unable to write to " + path;
    }

    sam_close(file);

    file = sam_open(path.c_str(), "r");

    int code = sam_index_build(path.c_str(), 0);
    if (code != 0) {
        throw "Cannot index " + path;
    }

    sam_close(file);
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}

    open_samFile_t(samFile* file, bam_hdr_t* header, hts_idx_t* idx) : file(file), header(header), idx(idx) {}
};

open_samFile_t* open_samFile(std::string fname_str, bool index_file = false) {
    const char* fname = fname_str.c_str();
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
	if (f) {
		hts_idx_destroy(f->idx);
		bam_hdr_destroy(f->header);
		sam_close(f->file);
		delete f;
	}
}

struct bam_redux_t {
    static const uint8_t IS_REV = 1, IS_MREV = 2;

    hts_pos_t start, end, mstart, isize;
    int left_clip_size, right_clip_size;
    uint8_t flag = 0;
    std::vector<uint8_t> seq;
    std::vector<uint8_t> qual;
    std::vector<uint32_t> cigar;

    bam_redux_t() {}
    bam_redux_t(bam1_t* read) : start(read->core.pos), end(bam_endpos(read)), mstart(read->core.mpos),
    isize(read->core.isize), left_clip_size(get_left_clip_size(read)), right_clip_size(get_right_clip_size(read)) {

        if (bam_is_rev(read)) flag |= IS_REV;
        if (bam_is_mrev(read)) flag |= IS_MREV;

        uint8_t* seq_array = bam_get_seq(read);
        seq = std::vector<uint8_t>(seq_array, seq_array+(read->core.l_qseq+1)/2);

        uint8_t* qual_array = bam_get_qual(read);
        qual = std::vector<uint8_t>(qual_array, qual_array+read->core.l_qseq);

        uint32_t* cigar_array = bam_get_cigar(read);
        cigar = std::vector<uint32_t>(cigar_array, cigar_array+read->core.n_cigar);
    }

    int seq_len() {
        return qual.size();
    }

    bool is_rev() {
        return flag & IS_REV;
    }
    bool is_mrev() {
        return flag & IS_MREV;
    }

    hts_pos_t unclipped_start() {
        return start-left_clip_size;
    }
    hts_pos_t unclipped_end() {
        return end+right_clip_size;
    }

    std::string cigar_string() {
        std::stringstream ss;
        for (uint32_t c : cigar) ss << bam_cigar_oplen(c) << bam_cigar_opchr(c);
        return ss.str();
    }
};

#endif //SMALLINSFINDER_SAM_UTILS_H
