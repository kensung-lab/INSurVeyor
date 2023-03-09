#ifndef VCF_UTILS_H_
#define VCF_UTILS_H_

#include <ctime>
#include <chrono>
#include <sstream>

#include "htslib/vcf.h"
#include "utils.h"

bcf_hrec_t* generate_contig_hrec() {
	bcf_hrec_t* contig_hrec = new bcf_hrec_t;
	contig_hrec->type = BCF_HL_CTG;
	contig_hrec->key = strdup("contig");
	contig_hrec->value = NULL;
	contig_hrec->keys = contig_hrec->vals = NULL;
	contig_hrec->nkeys = 0;
	int r1 = bcf_hrec_add_key(contig_hrec, "ID", 2);
	int r2 = bcf_hrec_add_key(contig_hrec, "length", 6);
	if (r1 || r2) {
		throw std::runtime_error("Failed to create contig to VCF header.");
	}
	return contig_hrec;
}
bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = bcf_hdr_init("w");

	// add contigs
	for (std::string contig_name : contigs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(contigs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(header, hrec);
	}

	int len;

	// add FILTER
	char size_flt_tag[1000];
	sprintf(size_flt_tag, "##FILTER=<ID=SMALL,Description=\"Insertion smaller than %d bp.\">", config.min_insertion_size);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, size_flt_tag, &len));

	const char* anomalous_sc_flt_tag = "##FILTER=<ID=ANOMALOUS_SC_NUMBER,Description=\"The number of soft-clipped reads supporting this call is too large.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_sc_flt_tag, &len));

	const char* anomalous_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_DEPTH,Description=\"The insertion region has anomalous depth.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_depth_flt_tag, &len));

	const char* alt_short_flt_tag = "##FILTER=<ID=ALT_SHORTER_THAN_REF,Description=\"If this insertion/replacement was real, alternative"
			"allele would be shorter than reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, alt_short_flt_tag, &len));

	const char* low_support_flt_tag = "##FILTER=<ID=LOW_SUPPORT,Description=\"Insertion has low support.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_support_flt_tag, &len));

	const char* no_disc_support_flt_tag = "##FILTER=<ID=NO_DISC_SUPPORT,Description=\"Insertion has no discordant support.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, no_disc_support_flt_tag, &len));

	const char* low_score_flt_tag = "##FILTER=<ID=LOW_SCORE,Description=\"Evidence against the insertion overwhelms the evidence in its favour.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_score_flt_tag, &len));

	const char* homopolymer_flt_tag = "##FILTER=<ID=HOMOPOLYMER_INSSEQ,Description=\"Inserted sequence is a homopolymer.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, homopolymer_flt_tag, &len));

	const char* mh_long_flt_tag = "##FILTER=<ID=MH_TOO_LONG,Description=\"Microhomology at the breakpoints is too long.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mh_long_flt_tag, &len));

	const char* lq_ass_flt_tag = "##FILTER=<ID=LOW_QUAL_ASSEMBLY,Description=\"Assembled quality is deemed not reliable.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, lq_ass_flt_tag, &len));

	const char* not_enoug_dp_flt_tag = "##FILTER=<ID=NOT_ENOUGH_DISC_PAIRS,Description=\"Not enough discordant pairs support this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enoug_dp_flt_tag, &len));

	// add INFO tags
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, end_tag, &len));

	const char* mhlen_tag = "##INFO=<ID=MHLEN,Number=1,Type=Integer,Description=\"Microhomology length.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, mhlen_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svlen_tag, &len));

	const char* svinslen_tag = "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinslen_tag, &len));

	const char* svinsseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinsseq_tag, &len));

	const char* overlap_tag = "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Overlap (in bp) between the left and right contigs.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, overlap_tag, &len));

	const char* disc_tag = "##INFO=<ID=DISCORDANT,Number=2,Type=Integer,Description=\"Discordant pairs supporting the left and right breakpoints of this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, disc_tag, &len));

	const char* sr_tag = "##INFO=<ID=SPLIT_READS,Number=2,Type=Integer,Description=\"Split reads supporting the left and right breakpoints of this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sr_tag, &len));

	const char* fwd_sr_tag = "##INFO=<ID=FWD_SPLIT_READS,Number=2,Type=Integer,Description=\"Forward split reads supporting the left and right breakpoints of this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, fwd_sr_tag, &len));

	const char* rev_sr_tag = "##INFO=<ID=REV_SPLIT_READS,Number=2,Type=Integer,Description=\"Reverse split reads supporting the left and right breakpoints of this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rev_sr_tag, &len));

	const char* span_tag = "##INFO=<ID=SPANNING_READS,Number=2,Type=Integer,Description=\"Negative evidence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, span_tag, &len));

	const char* scores_tag = "##INFO=<ID=SCORES,Number=2,Type=Float,Description=\"Scores for the the left and right breakpoints of this insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, scores_tag, &len));

	const char* left_anchor_tag = "##INFO=<ID=LEFT_ANCHOR,Number=1,Type=String,Description=\"Region flanking the insertion to the left, that contains"
			"reads supporting the insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, left_anchor_tag, &len));

	const char* lac_tag = "##INFO=<ID=LEFT_ANCHOR_CIGAR,Number=1,Type=String,Description=\"CIGAR of the alignment between the reference"
			"and LEFT ANCHOR.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, lac_tag, &len));

	const char* rac_tag = "##INFO=<ID=RIGHT_ANCHOR_CIGAR,Number=1,Type=String,Description=\"CIGAR of the alignment between the reference"
				"and RIGHT ANCHOR.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, rac_tag, &len));

	const char* right_anchor_tag = "##INFO=<ID=RIGHT_ANCHOR,Number=1,Type=String,Description=\"Region flanking the insertion to the right, that contains"
			"reads supporting the insertion.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, right_anchor_tag, &len));

	const char* tqc_tag = "##INFO=<ID=TRANS_QUERY_COV,Number=2,Type=Integer,Description=\"Length of the prefix and suffix of the transposed sequence that was covered by reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, tqc_tag, &len));

	const char* sd_tag = "##INFO=<ID=STABLE_DEPTHS,Number=2,Type=Integer,Description=\"Depths of the stable regions (in practice, the regions left and right of the insertion site).\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, sd_tag, &len));

	const char* avg_nm_tag = "##INFO=<ID=AVG_STABLE_NM,Number=2,Type=Integer,Description=\"Average NM of stable reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, avg_nm_tag, &len));

	const char* algo_tag = "##INFO=<ID=ALGORITHM,Number=1,Type=String,Description=\"Algorithm used to report the call.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, algo_tag, &len));

	const char* precise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Breakpoints and/or inserted sequence are imprecise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, precise_tag, &len));

	const char* incomplete_ass_tag = "##INFO=<ID=INCOMPLETE_ASSEMBLY,Number=0,Type=Flag,Description=\"Inserted sequence is too long and only part of it could be assembled.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, incomplete_ass_tag, &len));

	const char* full_junc_seq_tag = "##INFO=<ID=FULL_JUNCTION_SEQ,Number=0,Type=Flag,Description=\"Full junction sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, full_junc_seq_tag, &len));

	const char* src_reg_seq_tag = "##INFO=<ID=SOURCE_REGION,Number=0,Type=Flag,Description=\"Source region chosen by the one-end remapping.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, src_reg_seq_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, gt_tag, &len));

	// add ALT
	const char* ins_alt_tag = "##ALT=<ID=INS,Description=\"Insertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ins_alt_tag, &len));

	std::string cmd_tag = "##INSurVeyorCommand=" + command;
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##INSurVeyorVersion=" + config.version + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, version_tag.c_str(), &len));

	std::stringstream called_by_ss;
	called_by_ss << "##calledBy=INSurVeyor " << config.version << "; ";
	called_by_ss << "seed: " << config.seed << "; ";
	called_by_ss << "max-clipped-pos-dist: " << config.max_clipped_pos_dist << "; ";
	called_by_ss << "min-insertion-size: " << config.min_insertion_size << "; ";
	called_by_ss << "max-trans-size: " << config.max_trans_size << "; ";
	called_by_ss << "min-stable-mapq: " << config.min_stable_mapq << "; ";
	called_by_ss << "min-clip-len: " << config.min_clip_len << "; ";
	called_by_ss << "max-seq-error: " << config.max_seq_error << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	called_by_ss << "per-contig-stats: " << (config.per_contig_stats ? "true" : "false") << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	// add samples
	bcf_hdr_add_sample(header, sample_name.c_str());

	return header;
}

int get_sv_end(bcf1_t* sv, bcf_hdr_t* hdr) {
	bcf_unpack(sv, BCF_UN_INFO);

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

    throw std::runtime_error("SV " + std::string(sv->d.id) + " has no END or SVLEN annotation.");
}

std::string get_ins_seq(bcf1_t* sv, bcf_hdr_t* hdr) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	bcf_unpack(sv, BCF_UN_INFO);
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

std::string get_left_anchor(bcf1_t* sv, bcf_hdr_t* hdr) {
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "LEFT_ANCHOR", (void**) &data, &size);
	if (data) return data;
	return "";
}
std::string get_right_anchor(bcf1_t* sv, bcf_hdr_t* hdr) {
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "RIGHT_ANCHOR", (void**) &data, &size);
	if (data) return data;
	return "";
}

void insertion_to_bcf_entry(insertion_t* insertion, bcf_hdr_t* hdr, bcf1_t* bcf_entry, std::string id, chr_seqs_map_t& contigs) {
	bcf_clear(bcf_entry);
	bcf_entry->rid = bcf_hdr_name2id(hdr, insertion->chr.c_str());
	bcf_entry->pos = insertion->start;
	bcf_update_id(hdr, bcf_entry, id.c_str());

	std::string alleles = std::string(1, contigs.get_seq(insertion->chr)[insertion->start]) + ",<INS>";
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	int int_conv; // current htslib does not support writing int64 yet

	int_conv = insertion->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", "INS");
	if (insertion->ins_seq.find("-") == std::string::npos) {
		int_conv = insertion->ins_seq.length();
		bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		int_conv = insertion->ins_seq.length() - (insertion->end-insertion->start);
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
	} else {
		bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
	}
	bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", insertion->ins_seq.c_str());

	int int2_conv[2];
	int2_conv[0] = insertion->rc_reads(), int2_conv[1] = insertion->lc_reads();
	bcf_update_info_int32(hdr, bcf_entry, "SPLIT_READS", int2_conv, 2);
	int2_conv[0] = insertion->rc_fwd_reads, int2_conv[1] = insertion->lc_fwd_reads;
	bcf_update_info_int32(hdr, bcf_entry, "FWD_SPLIT_READS", int2_conv, 2);
	int2_conv[0] = insertion->rc_rev_reads, int2_conv[1] = insertion->lc_rev_reads;
	bcf_update_info_int32(hdr, bcf_entry, "REV_SPLIT_READS", int2_conv, 2);

	// add GT info
	int gt[1];
	gt[0] = bcf_gt_unphased(1);
	bcf_update_genotypes(hdr, bcf_entry, gt, 1);
}

#endif /* VCF_UTILS_H_ */
