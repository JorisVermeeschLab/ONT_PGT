#!/bin/bash -l
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="10"
#SBATCH --ntasks="10"
#SBATCH --mem-per-cpu=5G
#SBATCH --time="1-00:00:00"
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"

# sbatch --cluster=genius mapping_CHM13v1.0_hifi.slurm

SAMPLE="xx"
RG="HG002"
PROJECT_DIR="xx";
fastq_DIR="$PROJECT_DIR/xx_nanofilt_trimmed";
wkdir="$PROJECT_DIR/xx"
OUTPUT_DIR="$wkdir/1_mapped";

REF="/lustre1/project/stg_00019/research/yan/nanopore_data/00_06_T2T-CHM13v1.0/chm13.draft_v1.0.fasta"

# minimap2 can mapping PacBio or Oxford Nanopore genomic reads to the human genome;
# needs to be tuned; choose a preset with option -x
# --> default setting is the same as map-ont

module load minimap2/2.12-foss-2018a
module load SAMtools/1.9-GCC-6.4.0-2.28

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR

#  -a output in the SAM format
#  -R STR SAM read group line in a format like '@RG\tID:foo\tSM:bar' 

minimap2 -ax map-pb -R "@RG\tID:$RG\tSM:$RG\tPL:$RG\tLB:$RG\tPU:$RG" -t 20 $REF $fastq_DIR/xx.fastq.gz > $SAMPLE.sam

# Convert to bam file
samtools view -Shb $SAMPLE.sam > $SAMPLE.bam
# rm $SAMPLE.sam 

# Sort the bam file
samtools sort -o $SAMPLE.sorted.bam $SAMPLE.bam
rm $SAMPLE.bam

# create an index file
samtools index $SAMPLE.sorted.bam

