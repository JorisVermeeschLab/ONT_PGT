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
#SBATCH --time="24:00:00"

# submit job
# sbatch --cluster=genius GRCh38_rm_lowqual_phasing.slurm

module load tabix/0.2.6-GCCcore-6.4.0

WKDIR="xx"
cd $WKDIR
# exlude low quality callls
zcat merge_output.vcf.gz|grep '^#'> filtered.merge_output.vcf
zcat merge_output.vcf.gz|grep -v '#'|awk -F: '$(NF-2)>2'>> filtered.merge_output.vcf

# compress and creat index 

bgzip filtered.merge_output.vcf
tabix filtered.merge_output.vcf.gz

# phasing (single)
WHATSHAP_DIR="/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/bin"
OUT_DIR=$WKDIR/whatshap_phasing
BAM="xx"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"

mkdir -p $OUT_DIR

$WHATSHAP_DIR/whatshap phase \
--ignore-read-groups \
-o $OUT_DIR/phased.filtered.merge_output.vcf.gz \
--reference=$REF \
filtered.merge_output.vcf.gz $BAM

# creat index
tabix $OUT_DIR/phased.filtered.merge_output.vcf.gz

# get phasing stat & gtf file for phase set
$WHATSHAP_DIR/whatshap stats --gtf=$OUT_DIR/phased.gtf $OUT_DIR/phased.filtered.merge_output.vcf.gz>$OUT_DIR/phasing_statistics.txt
