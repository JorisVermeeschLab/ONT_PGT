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
#SBATCH --ntasks-per-node="5"
#SBATCH --ntasks="5"
#SBATCH --time="30:00:00"

# submit job
# sbatch --cluster=genius 6.2_phase_snp_indel.slurm

module load tabix/0.2.6-GCCcore-6.4.0

WKDIR="xx"
cd $WKDIR

# phasing (single)
WHATSHAP_DIR="/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/bin"
OUT_DIR=$WKDIR/whatshap_phasing_include_indel
BAM=""
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"

mkdir -p $OUT_DIR

$WHATSHAP_DIR/whatshap phase \
--ignore-read-groups \
-o $OUT_DIR/phased.filtered.merge_output.vcf.gz \
--indels \
--reference=$REF \
filtered.merge_output.vcf.gz $BAM

# creat index
tabix $OUT_DIR/phased.filtered.merge_output.vcf.gz

# get phasing stat & gtf file
$WHATSHAP_DIR/whatshap stats --gtf=$OUT_DIR/phased.gtf $OUT_DIR/phased.filtered.merge_output.vcf.gz>$OUT_DIR/phasing_statistics.txt
