#!/bin/bash
# chmod 755 1_merge_vcf_and_extract_required_info.sh
# nohup 1_merge_vcf_and_extract_required_info.sh>1_merge_vcf_and_extract_required_info.sh.log 2>&1 &
# r23i27n24
# 21686

module load BCFtools/1.9-foss-2018a

wkdir=/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_1_2_test_explore_PGT_M_merge_whole_vcf


# -------path to parental phased vcf files---------------------

# HG003_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_03_resource_parental_pacbio_data/HG003.phased.vcf.gz"
# HG004_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_03_resource_parental_pacbio_data/HG004.phased.vcf.gz"

# # -------path to son data phased vcf files---------------------

# bulk_no_10x_LSK110="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/0_bulk_no_10x_LSK110/6_variant_calling_guppy4_model"
# sc_PTA_10x_LSK110="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/2_sc_PTA_10x_LSK110/6_variant_calling_new_model"
# mc_PTA_4x_LSK110="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/4_mc_PTA_4x_LSK110/6_variant_calling_new_model"
# sc_MDA_24x_LSK110="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/7_sc_MDA_24x_LSK110/6_variant_calling_new_model"
# mc_MDA_23x_LSK110="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_data_LSK110_HG002/8_mc_MDA_23x_LSK110/6_variant_calling_new_model"
# sc_MDA_22x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/2_sc_MDA_22X_LSK114/6_new_model_HG002_SC_LSK114_variant_calling"
# mc_MDA_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/1_mc_MDA_24x_LSK114/6_new_model_HG002_MC_LSK114_variant_calling"
# bulk_no_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/3_bulk_no_24x_LSK114/6_variant_calling"

# -------merge vcf files---------------------

for data in "bulk_no_10x_LSK110" "sc_PTA_10x_LSK110" "mc_PTA_4x_LSK110" "sc_MDA_24x_LSK110" "mc_MDA_23x_LSK110" "sc_MDA_22x_LSK114" "mc_MDA_24x_LSK114" "bulk_no_24x_LSK114"
do
cd $wkdir/$data
# link=`echo "${!data}"`

# # merge 3 vcf files 
# bcftools merge $link/phased_merge_output.vcf.gz  $HG003_phased_vcf  $HG004_phased_vcf > "$data"_merged.vcf 

# # extract  info from merged vcf 
# echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG002"\t"GQ_HG003"\t"GQ_HG004"\t"DP_HG002"\t"DP_HG003"\t"DP_HG004"\t"PS_HG002"\t"PS_HG003"\t"PS_HG004 > "$data"_merged_vcf_info.txt
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' "$data"_merged.vcf >> "$data"_merged_vcf_info.txt

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG003&HG004
# 3. 4. remove region where HG003==HG004==hom
# 5. 6. replace 0/0 1/1 with 0|0 1|1 
# 7.8.9.  replace 0/1 1/0 ./. in HG002 with .|.

head -1 "$data"_merged_vcf_info.txt>"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|awk '$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$6!="0/1"&&$7!="0/1"' \
|awk '!($6=="1/1"&&$7=="1/1")'\
|awk '!($6=="0/0"&&$7=="0/0")'\
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'\
|sed 's#0/1#.|.#g'\
|sed 's#1/0#.|.#g'\
|sed 's#./.#.|.#g'>> "$data"_merged_vcf_info_filtered_revised.txt

done