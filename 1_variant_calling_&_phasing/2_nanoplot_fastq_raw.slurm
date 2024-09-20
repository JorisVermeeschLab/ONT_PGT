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
#SBATCH --ntasks-per-node="6"
#SBATCH --ntasks="6"
#SBATCH --time="7:00:00"

# chaek quality of the merged raw fastq file with NanoPlot

# submit job
# sbatch --cluster=genius 2_nanoplot_fastq_raw.slurm

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load Python/3.7.2-foss-2018a  

SAMPLE=""
PROJECT_DIR="";
SAMPLE_DIR="$PROJECT_DIR/1_merged_fastq";
OUTPUT_DIR="$PROJECT_DIR/2_nanoplot_fastq_raw";
mkdir -p $OUTPUT_DIR;

conda activate nanoplot

NanoPlot --fastq $SAMPLE_DIR/$SAMPLE.fastq.gz \
    -o $OUTPUT_DIR \
    -t 18