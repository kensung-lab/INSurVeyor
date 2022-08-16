# INSurVeyor
An insertion caller for Illumina paired-end WGS data.

## Installation

There are two options for obtaining. If singularity is available on your system, you can download the image in Releases. Alternatively, you can compile the code yourself. This should be 
straightforward, but some dependencies must be installed.

### Singularity

The latest image can be obtained under Releases.

### Compiling the code

The following commands should be sufficient:

```
git clone https://github.com/kensung-lab/INSurVeyor
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

If htslib does not build correctly, please refer to https://github.com/samtools/htslib

Python is necessary to run INSurVeyor. Libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. If 
Python 2 is used, numpy 1.16.6, pyfaidx 0.5.9 and pysam 0.16.0.1 are the recommended (i.e., tested) versions. If Python 3 is used, then numpy 1.21.2, pyfaidx 0.5.9.1 and pysam 0.16.0.1 were 
tested.

## Running

INSurVeyor needs a BAM/CRAM file, a (possibly empty) working directory and reference genome in FASTA format.
The BAM/CRAM file must be coordinate-sorted and indexed. Furthermore, the MC and the MQ tag must be present for all primary alignments, when applicable.

Recent versions of BWA MEM (0.7.17) will add the MC tag. The easiest (but probably not the fastest) way to add the MQ tag is to use Picard FixMateInformation 
(http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation) 
```
java -jar picard.jar FixMateInformation I=file.bam
```

If you are using the singularity image, you can run INSurVeyor with the following command:
```
singularity run insurveyor.sif --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA
```

The command
```
singularity run insurveyor.sif
```
will print the help message of the software.

If you compiled the code, you can use 
```
python surveyor.py --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA
```

## Output

The output is a standard VCF file. It will be placed under WORKDIR/out.pass.vcf.gz. These are the insertions that INSurVeyor deemed confident enough. 

The file WORKDIR/out.vcf.gz contains all of the insertions, including those that did not pass the filters. Most of them will be false positives. It is not recommend to you use this file unless 
for specific situations (e.g., you are looking for something specific).

## Citation

A manuscript on INSurVeyor is currently in preparation.

