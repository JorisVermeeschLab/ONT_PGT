#!/bin/bash
# chmod 755 1_merge_vcf_and_extract_required_info.sh
# nohup 1_merge_vcf_and_extract_required_info.sh>1_merge_vcf_and_extract_required_info.sh.log 2>&1 &
# r23i27n22
# 

module load BCFtools/1.9-foss-2018a

# GRCh38 1_check_GW_PGT_M

## ------------------------merge vcf files-------------------------
data_dir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/GRCh38/1_check_GW_PGT_M"
category_1="independent_vc_and_independent_phasing"
category_2="trio_vc_and_pedgree_phasing"

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/GRCh38/2_check_multiple_filtering"

## category_1 independent_vc_and_independent_phasing

for data in "independent_sc_PTA_10x_LSK110" "independent_mc_PTA_4x_LSK110" "independent_sc_MDA_24x_LSK110" "independent_mc_MDA_23x_LSK110" "independent_mc_MDA_24x_LSK114" "independent_sc_MDA_22x_LSK114" "independent_bulk_no_24x_LSK114" 
do

cd $data_dir/$category_1/$data

outdir=$wkdir/$category_1/$data
mkdir -p $outdir

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG002&HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG002&HG003&HG004
# 3. 4. replace 0/0 1/1 with 0|0 1|1 

head -1 "$data"_merged_vcf_info.txt> $outdir/"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|awk '$5!="./."&&$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$5!="0/1"&&$6!="0/1"&&$7!="0/1"' \
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> $outdir/"$data"_merged_vcf_info_filtered_revised.txt;

done

## category_2 trio_vc_and_pedgree_phasing

for data in "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_23x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do

cd $data_dir/$category_2/$data

outdir=$wkdir/$category_2/$data
mkdir -p $outdir

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG002&HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG002&HG003&HG004
# 3. 4. replace 0/0 1/1 with 0|0 1|1 

head -1 "$data"_merged_vcf_info.txt> $outdir/"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|awk '$5!="./."&&$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$5!="0/1"&&$6!="0/1"&&$7!="0/1"' \
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> $outdir/"$data"_merged_vcf_info_filtered_revised.txt;

done