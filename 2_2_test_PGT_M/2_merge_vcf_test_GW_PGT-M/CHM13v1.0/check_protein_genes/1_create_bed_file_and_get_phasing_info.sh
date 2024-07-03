#!/bin/bash

# chmod 755 1_create_bed_file_and_get_phasing_info.sh
# nohup 1_create_bed_file_and_get_phasing_info.sh>1_create_bed_file_and_get_phasing_info.sh.log 2>&1 &
# r23i27n22
# 32980

peddir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/1_check_GW_PGT_M/trio_vc_and_pedgree_phasing"

uniq_gene_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/T2T-CHM13v1.0/noYM_uniq_gene_protein_coding.sorted.merged.bed"

module load BEDTools
module load BCFtools

for data in "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do

vcf=$peddir/$data/"$data"_merged.vcf

outdir="$peddir/check_protein_genes/$data"
mkdir -p $outdir
cd $outdir

bedtools intersect -header -a $vcf \
                             -b $uniq_gene_protein_coding_bed> $data.protein_genes.vcf

# extract  info from merged vcf 
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG002"\t"GQ_HG003"\t"GQ_HG004"\t"DP_HG002"\t"DP_HG003"\t"DP_HG004"\t"PS_HG002"\t"PS_HG003"\t"PS_HG004 > $data.protein_genes.vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' $data.protein_genes.vcf >>  $data.protein_genes.vcf_info.txt

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG002&HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG002&HG003&HG004
# 3. 4. remove region where HG003==HG004==hom
# 5. 6. replace 0/0 1/1 with 0|0 1|1 

head -1 $data.protein_genes.vcf_info.txt>$data.protein_genes.vcf_filtered_revised.txt

cat $data.protein_genes.vcf_info.txt\
|awk '$5!="./."&&$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$5!="0/1"&&$6!="0/1"&&$7!="0/1"' \
|awk '!($6=="1/1"&&$7=="1/1")'\
|awk '!($6=="0/0"&&$7=="0/0")'\
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> $data.protein_genes.vcf_filtered_revised.txt;

done