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
#SBATCH --ntasks-per-node="10"
#SBATCH --ntasks="10"
#SBATCH --time="10:00:00"

# submit job
# sbatch --cluster=genius 4_nanoplot_fastq_trimmed.slurm

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load Python/3.7.2-foss-2018a  

SAMPLE=""
PROJECT_DIR="";
SAMPLE_DIR="$PROJECT_DIR/3_nanofilt_trimmed";
OUTPUT_DIR="$PROJECT_DIR/4_nanoplot_fastq_trimmed";
mkdir -p $OUTPUT_DIR;

conda activate nanoplot

fastq=`ls $SAMPLE_DIR/*.fastq.gz`

NanoPlot --fastq $fastq \
    -o $OUTPUT_DIR \
    -t 18