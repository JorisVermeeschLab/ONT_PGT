#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --get-user-env="60L"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="12"
#SBATCH --ntasks="12"
#SBATCH --time="30:00:00"

# submit job
# sbatch --cluster=genius 5_mapping.slurm

SAMPLE=""
PROJECT_DIR="";

SAMPLE_DIR="$PROJECT_DIR/3_nanofilt_trimmed";
OUTPUT_DIR="$PROJECT_DIR/5_mapped";
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"


module load minimap2/2.12-foss-2018a
module load SAMtools/1.9-GCC-6.4.0-2.28

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR

#  -a output in the SAM format
#  -R STR SAM read group line in a format like '@RG\tID:foo\tSM:bar' 

minimap2 -a -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" -t 20 $REF $SAMPLE_DIR/$SAMPLE.trimmed.fastq.gz | samtools sort > $SAMPLE.trimmed.sorted.bam

# create an index file
samtools index $SAMPLE.trimmed.sorted.bam

