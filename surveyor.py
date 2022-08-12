from __future__ import print_function
import argparse, os, glob, time
import pysam, pyfaidx
from random_pos_generator import RandomPositionGenerator
import numpy as np

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

cmd_parser = argparse.ArgumentParser(description='INSurVeyor, an insertion caller.')
cmd_parser.add_argument('bam_file', help='Input bam file.')
cmd_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--max_clipped_pos_dist', type=int, default=5, help='Max distance (in bp) for two clips to be considered '
                                                                   'representing the same breakpoint.')
cmd_parser.add_argument('--min_insertion_size', type=int, default=50, help='Minimum size of the insertion to be called.'
                                                                     'Smaller insertions will be reported anyway but marked'
                                                                     'as SMALL in the FILTER field.')
cmd_parser.add_argument('--max_insertion_size', type=int, default=10000, help='Maximum size of the insertions which '
                                                                              'INSurVeyor will try to predict.')
cmd_parser.add_argument('--min_stable_mapq', type=int, default=20, help='Minimum MAPQ for a stable read.')
cmd_parser.add_argument('--min_clip_len', type=int, default=15, help='Minimum clip len to consider.')
cmd_parser.add_argument('--max_seq_error', type=float, default=0.04, help='Max sequencing error admissible on the platform used.')
cmd_parser.add_argument('--sampling-regions', help='File in BED format containing a list of regions to be used to estimate'
                                                   'statistics such as depth.')
cmd_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                                 'If not provided, the basename of the bam/cram file will be used,'
                                                                 'up until the first \'.\'')
cmd_parser.add_argument('--sample-clipped-pairs', action='store_true', help='When estimating the insert size distribution '
                                                                            'by sampling pairs, do not discard pairs where '
                                                                            'one or both of the reads are clipped.')
cmd_args = cmd_parser.parse_args()

SURV_PATH = os.path.dirname(os.path.realpath(__file__))

# Create config file in workdir
config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("max_clipped_pos_dist %d\n" % cmd_args.max_clipped_pos_dist)
config_file.write("min_insertion_size %d\n" % cmd_args.min_insertion_size)
config_file.write("max_insertion_size %d\n" % cmd_args.max_insertion_size)
config_file.write("min_stable_mapq %d\n" % cmd_args.min_stable_mapq)
config_file.write("min_clip_len %d\n" % cmd_args.min_clip_len)
config_file.write("max_seq_error %f\n" % cmd_args.max_seq_error)

# Find read length
read_len = 0
bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
for i, read in enumerate(bam_file.fetch(until_eof=True)):
    if i > MAX_READS: break
    read_len = max(read_len, read.query_length)
config_file.write("read_len %d\n" % read_len)


contig_map = open("%s/contig_map" % cmd_args.workdir, "w")
for i, k in enumerate(bam_file.references):
    contig_map.write("%s\n" % (k));
contig_map.close();

# Generate general distribution of insert sizes
reference_fa = pyfaidx.Fasta(cmd_args.reference)
rand_pos_gen = RandomPositionGenerator(reference_fa, cmd_args.sampling_regions)
random_positions = []
for i in range(1,1000001):
    if i % 100000 == 0: print(i, "random positions generated.")
    random_positions.append(rand_pos_gen.next())

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)

general_dist = []
avg_depth = 0
samplings = 0
rnd_i = 0
while rnd_i < len(random_positions) and len(general_dist) < GEN_DIST_SIZE:
    chr, pos = random_positions[rnd_i]
    rnd_i += 1

    if pos > len(reference_fa[chr])-10000:
        continue

    samplings += 1
    i = 0
    for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
        if read.reference_start < pos:
            avg_depth += 1
        if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
        0 < read.template_length < MAX_ACCEPTABLE_IS and \
        (cmd_args.sample_clipped_pairs or 'S' not in read.cigarstring and 'S' not in read.get_tag('MC')):
            if i > 100: break
            i += 1
            general_dist.append(read.template_length)
reference_fa.close()

mean_is = np.mean(general_dist)
stddev_is = np.std(general_dist)
avg_depth = float(avg_depth)/samplings

print("Average depth:", avg_depth)

general_dist = [x for x in general_dist if abs(x-mean_is) < 5*stddev_is]

mean_is = int(np.mean(general_dist))
lower_stddev_is = int(np.sqrt(np.mean([(mean_is-x)**2 for x in general_dist if x < mean_is])))
higher_stddev_is = int(np.sqrt(np.mean([(x-mean_is)**2 for x in general_dist if x > mean_is])))

# TODO should move this to stats
min_is, max_is = mean_is-3*lower_stddev_is, mean_is+3.5*higher_stddev_is
config_file.write("min_is %d\n" % min_is)
config_file.write("avg_is %d\n" % mean_is)
config_file.write("max_is %d\n" % max_is)
config_file.write("read_len %d\n" % read_len)
config_file.close();

workspace = cmd_args.workdir + "/workspace"
if not os.path.exists(workspace):
    os.makedirs(workspace)
    os.makedirs(workspace + "/clipped/")
    os.makedirs(workspace + "/clip_consensuses/")
    os.makedirs(workspace + "/mateseqs/")
    os.makedirs(workspace + "/R/")
    os.makedirs(workspace + "/L/")

start = time.time()
read_categorizer_cmd = SURV_PATH + "/reads_categorizer %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
print("Executing:", read_categorizer_cmd)
os.system(read_categorizer_cmd)
end = time.time()
print("Reads categorized in %d [s]" % (end-start))

# for f in glob.glob(workspace + "/clipped/*.bam"):
#     prefix = f[:-4]
#     pysam.sort("-@", str(cmd_args.threads), "-o", "%s.cs.bam" % prefix, f)
#     os.rename("%s.cs.bam" % prefix, f)

start = time.time()
clip_consensus_builder_cmd = SURV_PATH + "/clip_consensus_builder %s" % (cmd_args.workdir)
print("Executing:", clip_consensus_builder_cmd)
os.system(clip_consensus_builder_cmd)
end = time.time()
print("Clip consensus built in %d [s]" % (end-start))

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

start = time.time()
call_insertions_cmd = SURV_PATH + "/call_insertions %s %s %s" % (cmd_args.workdir, cmd_args.reference, sample_name)
print("Executing:", call_insertions_cmd)
os.system(call_insertions_cmd)
end = time.time()
print("Small insertions called in %d [s]" % (end-start))

start = time.time()
dc_remapper_cmd = SURV_PATH + "/dc_remapper %s %s %s" % (cmd_args.workdir, cmd_args.reference, sample_name)
print("Executing:", dc_remapper_cmd)
os.system(dc_remapper_cmd)
end = time.time()
print("DC remapped in %d [s]" % (end-start))

start = time.time()
add_filtering_info_cmd = SURV_PATH + "/add_filtering_info %s %s %s" % (cmd_args.bam_file, cmd_args.workdir, cmd_args.reference)
print("Executing:", add_filtering_info_cmd)
os.system(add_filtering_info_cmd)
end = time.time()
print("Filtering info added in %d [s]" % (end-start))

filter_cmd = SURV_PATH + "/filter %s %s 0.25" % (cmd_args.workdir, cmd_args.reference)
os.system(filter_cmd)
