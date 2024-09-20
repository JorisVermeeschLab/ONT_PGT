#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --get-user-env="60L"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user=""
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="1"
#SBATCH --ntasks="1"
#SBATCH --time="7:00:00"

# check quality of merged raw fastq file

# submit job
# sbatch --cluster=genius 2_fastqc_old.slurm

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load FastQC/0.11.8-Java-1.8.0_162

#setting all parameters (these could be changed)
SAMPLE=""
PROJECT_DIR="";
SAMPLE_DIR="$PROJECT_DIR/1_merged_fastq";
OUTPUT_DIR="$PROJECT_DIR/2_fastqc_raw";
mkdir -p $OUTPUT_DIR;

cd $SAMPLE_DIR;

fastqc -o $OUTPUT_DIR -t 15 $SAMPLE.fastq.gz;

