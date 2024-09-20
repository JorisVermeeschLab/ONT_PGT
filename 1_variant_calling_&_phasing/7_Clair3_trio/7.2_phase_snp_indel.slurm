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
# sbatch --cluster=genius 7.2_phase_snp_indel.slurm

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BCFtools/1.9-foss-2018a
module load tabix/0.2.6-GCCcore-6.4.0

WHATSHAP_DIR="/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/envs/clair3/bin/"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"
SAMPLE=""
F=""
M=""

bam_sample="xx/$SAMPLE.trimmed.sorted.bam"
bam_f="xx/$F.trimmed.sorted.bam"
bam_m="xx/$M.trimmed.sorted.bam"

wkdir="xx/7_Clair3_trio"
cd $wkdir

outdir="$wkdir/whatshap_phasing_include_indel"
mkdir -p $outdir

# # -----step2. merge parental and biopsy's vcf files (nocall to ref)-----
filtered_F_vcf="xx/6_Clair3/filtered.merge_output.vcf.gz"
filtered_M_vcf="xx/6_Clair3/filtered.merge_output.vcf.gz"

bcftools merge --missing-to-ref filtered.$SAMPLE.vcf.gz $filtered_F_vcf $filtered_M_vcf> $outdir/pedgree_merged.vcf

# -------------------step3.pedgree phasing--------------------------- 
cd $outdir
bgzip pedgree_merged.vcf
tabix pedgree_merged.vcf.gz 

# create ped file using sample names in bam file 
echo -e ."\t"$SAMPLE"\t"$F"\t"$M"\t"."\t".> pedigree.ped

# pedgree phasing including indel
$WHATSHAP_DIR/whatshap phase --indels --ped pedigree.ped --reference=$REF -o pedgree.phased.vcf.gz pedgree_merged.vcf.gz $bam_sample $bam_f $bam_f
tabix pedgree.phased.vcf.gz 

rm pedgree_merged.vcf.gz*
#-----------------step2. get phased vcf for HG002 from pedgree.phased.vcf.gz---------------
bcftools view -Oz -s $SAMPLE pedgree.phased.vcf.gz > $SAMPLE.phased.vcf.gz
tabix $SAMPLE.phased.vcf.gz

$WHATSHAP_DIR/whatshap stats --gtf=$SAMPLE.phased.gtf $SAMPLE.phased.vcf.gz>$SAMPLE.phasing_statistics.txt