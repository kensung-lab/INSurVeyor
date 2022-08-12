#include <iostream>
#include <string>
#include <algorithm>

#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "utils.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;

std::string get_sv_type(bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_ins_seq(bcf1_t* sv) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		return sv->d.allele[1];
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) return data;

	return "";
}

bool normalise_ins(bcf1_t* vcf_record) {

	char* chr_seq = chr_seqs.get_seq(bcf_hdr_id2name(hdr, vcf_record->rid));
	int end = get_sv_end(vcf_record);

	// try to shorten deletion if ins_seq
	std::string ins_seq = get_ins_seq(vcf_record);
	std::string orig_ins_seq = ins_seq;
	std::replace(ins_seq.begin(), ins_seq.end(), ' ', '-');

	int start_is = 0;
	while (start_is < ins_seq.length() && vcf_record->pos < end && toupper(ins_seq[start_is]) == toupper(chr_seq[vcf_record->pos+1])) {
		start_is++;
		vcf_record->pos++;
	}
	ins_seq = ins_seq.substr(start_is);

	int end_is = ins_seq.length();
	while (end_is > 0 && vcf_record->pos < end && toupper(ins_seq[end_is-1]) == toupper(chr_seq[end])) {
		end_is--;
		end--;
	}
	ins_seq = ins_seq.substr(0, end_is);

	if (vcf_record->pos == end) {
		while (toupper(chr_seq[vcf_record->pos]) == toupper(ins_seq[ins_seq.length()-1])) {
			for (int i = ins_seq.length()-1; i >= 1; i--) {
				ins_seq[i] = ins_seq[i-1];
			}
			ins_seq[0] = toupper(chr_seq[vcf_record->pos]);
			vcf_record->pos--;
			end--;
		}
	}

	if (ins_seq != orig_ins_seq) {
		if (ins_seq.empty()) {
			bcf_update_info_string(hdr, vcf_record, "SVINSSEQ", NULL);
		} else {
			bcf_update_info_string(hdr, vcf_record, "SVINSSEQ", ins_seq.c_str());
			if (bcf_get_info_flag(hdr, vcf_record, "INCOMPLETE_ASSEMBLY", NULL, NULL)) {
				int len = ins_seq.length();
				bcf_update_info_int32(hdr, vcf_record, "SVLEN", &len, 1);
			}
		}
	}

	int end_1based = end+1;
	bcf_update_info_int32(hdr, vcf_record, "END", &end_1based, 1);

	std::string alleles = std::string(1, chr_seq[vcf_record->pos]) + "," + vcf_record->d.allele[1];
	bcf_update_alleles_str(hdr, vcf_record, alleles.c_str());

	return (vcf_record->pos >= 0 && ins_seq[ins_seq.length()-1] != '-');
}


bool normalise(bcf1_t* vcf_record) {
	std::string svtype = get_sv_type(vcf_record);
	if (svtype == "INS") {
		return normalise_ins(vcf_record);
	}
	return false;
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];

	chr_seqs.read_fasta_into_map(reference_fname);

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* vcf_record = bcf_init();

	std::vector<bcf1_t*> normalised_vcf_records;
	while (bcf_read(in_vcf_file, hdr, vcf_record) == 0) {
		bcf1_t* vcf_record_norm = bcf_dup(vcf_record);
		if (normalise(vcf_record_norm)) {
			normalised_vcf_records.push_back(vcf_record_norm);
		}
	}

	std::sort(normalised_vcf_records.begin(), normalised_vcf_records.end(),
			[](const bcf1_t* b1, const bcf1_t* b2) { return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos); });

	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write VCF header to " +  out_vcf_fname);
	}
	for (bcf1_t* vcf_record_norm : normalised_vcf_records) {
		if (bcf_write(out_vcf_file, hdr, vcf_record_norm) != 0) {
			throw std::runtime_error("Failed to write VCF record to " +  out_vcf_fname);
		}
	}

	hts_close(in_vcf_file);
	hts_close(out_vcf_file);

	tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
