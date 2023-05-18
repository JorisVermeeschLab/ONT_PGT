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
# sbatch --cluster=genius GRCh38_rm_lowqual_calls_add_ref_calls_pedgree_phasing.slurm

module load BEDTools/2.27.1-intel-2018a
module load VCFtools/0.1.16-foss-2018a-Perl-5.26.1
module load tabix/0.2.6-GCCcore-6.4.0
module load BCFtools/1.9-foss-2018a
module load SAMtools/1.9-GCC-6.4.0-2.28

wkdir="xx"
outdir="$wkdir/2_vcf_GQ_filtered_whatshap_pedgree_phasing_with_refcall"
mkdir -p $outdir
cd $outdir

# -----step1. add ref call from GQ filtered single variant calling to trio variant calling results-----

vcf_trio_HG002="$wkdir/xx"
filtered_vcf_with_ref="xx/filtered.merge_output.vcf.gz"

# 1. filter trio variant calling results with GQ>2 
zcat $vcf_trio_HG002|grep '^#'> filtered.trio.HG002.vcf
zcat $vcf_trio_HG002|grep -v '#'|awk -F: '$(NF-2)>2'>> filtered.trio.HG002.vcf

# 2. get ref calls from GQ filtered single calling 
zcat $filtered_vcf_with_ref|grep '^#'> ref.HG002.vcf
zcat $filtered_vcf_with_ref|grep '0/0'>> ref.HG002.vcf

# 3. output ref calls that do not overlap with filtered.trio.HG002.vcf
# -v Only report those entries in A that have no overlap in B.
bedtools intersect -v -a ref.HG002.vcf \
                             -b filtered.trio.HG002.vcf>ref.HG002.to_merge.vcf 

# 4. merge these ref calls with GQ filtered.trio.HG002.vcf 
cat ref.HG002.to_merge.vcf>>filtered.trio.HG002.vcf
vcf-sort --chromosomal-order filtered.trio.HG002.vcf|bgzip>filtered.ref_added.trio.HG002.vcf.gz 
tabix filtered.ref_added.trio.HG002.vcf.gz
# remove intermediate files
rm filtered.trio.HG002.vcf ref.HG002.vcf ref.HG002.to_merge.vcf

# -----step2. merge parental and son's vcf files-----

#AJ Father GM24149	HG003
#AJ Mother GM24143	HG004

filtered_HG003_vcf_with_ref="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24149/6_GRCh38/Clair3_refcall/filtered.merge_output.vcf.gz"
filtered_HG004_vcf_with_ref="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24143/6_GRCh38/2_Clair3_refcall/filtered.merge_output.vcf.gz"

# ALT !=. HG002!=./.:.:.:.
bcftools merge filtered.ref_added.trio.HG002.vcf.gz $filtered_HG003_vcf_with_ref $filtered_HG004_vcf_with_ref \
|awk -F "\t" '$5!="."&&$10!="./.:.:.:."{print}'>pedgree_merged_with_ref.vcf

WHATSHAP_DIR="/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/envs/clair3/bin/"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"
chr_len="/lustre1/project/stg_00019/research/yan/nanopore_data/hg38_contig_length.txt"
BAM_SON="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/7_sc_MDA_24x_LSK110/6_change_bamRG/sc_MDA_24x_LSK110.sorted.mod.rg.bam"
BAM_F="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24149/6_GRCh38/6_mapped/Father_GM24149_HG003.rg.changed.bam"
BAM_M="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24143/6_GRCh38/1_mapped/Mother_GM24143_HG004.rg.changed.bam"

name_HG002=`samtools view -H $BAM_SON|grep ^@RG|cut -f 5|awk -F "SM:" '{print $2}'`
name_HG003=`samtools view -H $BAM_F|grep ^@RG|cut -f 5|awk -F "SM:" '{print $2}'`
name_HG004=`samtools view -H $BAM_M|grep ^@RG|cut -f 5|awk -F "SM:" '{print $2}'`

echo -e $name_HG002"\n"$name_HG003"\n"$name_HG004>names.txt
bcftools reheader --samples names.txt pedgree_merged_with_ref.vcf> pedgree_merged_with_ref.reheader.vcf

# -------------------step2.pedgree phasing--------------------------- 
# create ped file using sample names in bam file 
echo -e ."\t"$name_HG002"\t"$name_HG003"\t"$name_HG004"\t"."\t".>pedigree.ped

# pedgree phasing
$WHATSHAP_DIR/whatshap phase --ped pedigree.ped --reference=$REF -o pedgree.phased.vcf.gz pedgree_merged_with_ref.reheader.vcf $BAM_SON $BAM_F $BAM_M
tabix pedgree.phased.vcf.gz 
rm pedgree_merged_with_ref.vcf pedgree_merged_with_ref.reheader.vcf

#-----------------step2. get phased vcf for HG002 from pedgree.phased.vcf.gz---------------
bcftools view -Oz -s $name_HG002 pedgree.phased.vcf.gz > HG002.phased.vcf.gz
tabix HG002.phased.vcf.gz

$WHATSHAP_DIR/whatshap stats --gtf=HG002.phased.gtf HG002.phased.vcf.gz>HG002.phasing_statistics.txt