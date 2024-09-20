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
#SBATCH --ntasks-per-node="1"
#SBATCH --ntasks="1"
#SBATCH --time="15:00:00"

# submit job
# sbatch --cluster=genius 3_nanofilt_trim.slurm

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load Python/3.7.2-foss-2018a  

SAMPLE=""
PROJECT_DIR="xx";
SAMPLE_DIR="$PROJECT_DIR/1_merged_fastq";
OUTPUT_DIR="$PROJECT_DIR/3_nanofilt_trimmed";

mkdir -p $OUTPUT_DIR;
cd $OUTPUT_DIR;
rsync -ahrL $SAMPLE_DIR/$SAMPLE.fastq.gz .;
gunzip *gz;

# remove all sequences shorter than 500 nucleotides (option -l)
# trim the first 9 nucleotides off all reads (option –headcrop)
# -q Filter on a minimum average read quality score
# SNP identification: you don't want too many low quality bases

/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/bin/NanoFilt -q 9 -l 500 --headcrop 9 < $SAMPLE.fastq|gzip > $SAMPLE.trimmed.fastq.gz