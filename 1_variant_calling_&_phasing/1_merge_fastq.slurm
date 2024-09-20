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

# merge raw fastq files into a single fastq file

# submit job
# sbatch --cluster=genius 1_merge_fastq.slurm

sample=""
small_fastq="xx/fastq_pass"
workdir=""
outdir=$workdir/1_merged_fastq
mkdir -p $outdir

cat $small_fastq/*.fastq.gz > $outdir/$sample.fastq.gz