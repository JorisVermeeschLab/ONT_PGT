#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --job-name="2_bam_coverage_samtools.slurm"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="2"
#SBATCH --ntasks="2"
#SBATCH --time="10:00:00"

# sbatch --cluster=genius 2_bam_coverage_samtools.slurm


module load SAMtools/1.9-GCC-6.4.0-2.28

SAMPLE="xx"
bam_path="xx"

cd $bam_path

samtools depth -aa $SAMPLE.trimmed.sorted.bam > $SAMPLE.samtools_depth.txt
