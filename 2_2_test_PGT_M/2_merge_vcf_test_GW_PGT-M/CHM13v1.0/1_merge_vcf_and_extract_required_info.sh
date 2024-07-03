#!/bin/bash
# chmod 755 1_merge_vcf_and_extract_required_info.sh
# nohup 1_merge_vcf_and_extract_required_info.sh>1_merge_vcf_and_extract_required_info.sh.log 2>&1 &
# r23i27n22
# 2496

module load BCFtools/1.9-foss-2018a

# CHM13v1.0 1_check_GW_PGT_M

#-------path to phased parental vcf files (independent_vc_and_independent_phasing)---------------------

# AJ Father	GM24149	HG003
# AJ Mother	GM24143	HG004

HG003_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24149/6_CHM13v1.0/Clair3_refcall/whatshap_phasing/phased.filtered.merge_output.vcf.gz"
HG004_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24143/6_CHM13v1.0/Clair3_refcall/whatshap_phasing/phased.filtered.merge_output.vcf.gz"

## -------path to son data phased vcf files---------------------

# independent_vc_and_independent_phasing
independent_bulk_no_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/bulk_no_24x_LSK114/2_Clair3_refcall/whatshap_phasing/phased.filtered.merge_output.vcf.gz"
independent_mc_MDA_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/mc_MDA_24x_LSK114/2_Clair3_refcall/whatshap_phasing/phased.filtered.merge_output.vcf.gz"
independent_sc_MDA_22x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/sc_MDA_22X_LSK114/2_Clair3_refcall/whatshap_phasing/phased.filtered.merge_output.vcf.gz"

# trio_vc_and_pedgree_phasing
ped_bulk_no_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/bulk_no_24x_LSK114/2_Clair3_trio/vcf_GQ_filtered_whatshap_pedgree_phasing_with_refcall/HG002.phased.vcf.gz"
ped_mc_MDA_24x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/mc_MDA_24x_LSK114/2_Clair3_trio/vcf_GQ_filtered_whatshap_pedgree_phasing_with_refcall/HG002.phased.vcf.gz"
ped_sc_MDA_22x_LSK114="/lustre1/project/stg_00019/research/yan/nanopore_data/01_02_data_LSK114_HG002/map_to_CHM13v1.0/sc_MDA_22X_LSK114/2_Clair3_trio/vcf_GQ_filtered_whatshap_pedgree_phasing_with_refcall/HG002.phased.vcf.gz"

## ------------------------merge vcf files-------------------------
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/1_check_GW_PGT_M"
category_1="independent_vc_and_independent_phasing"
category_2="trio_vc_and_pedgree_phasing"

## category_1 independent_vc_and_independent_phasing

for data in "independent_mc_MDA_24x_LSK114" "independent_sc_MDA_22x_LSK114" "independent_bulk_no_24x_LSK114" 
do
HG002_phased_vcf=`echo "${!data}"`
mkdir $wkdir/$category_1/$data
cd $wkdir/$category_1/$data
# merge 3 vcf files 
bcftools merge $HG002_phased_vcf $HG003_phased_vcf $HG004_phased_vcf > "$data"_merged.vcf 

# extract  info from merged vcf 
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG002"\t"GQ_HG003"\t"GQ_HG004"\t"DP_HG002"\t"DP_HG003"\t"DP_HG004"\t"PS_HG002"\t"PS_HG003"\t"PS_HG004 > "$data"_merged_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' "$data"_merged.vcf >> "$data"_merged_vcf_info.txt

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG002&HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG002&HG003&HG004
# 3. 4. remove region where HG003==HG004==hom
# 5. 6. replace 0/0 1/1 with 0|0 1|1 

head -1 "$data"_merged_vcf_info.txt>"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|awk '$5!="./."&&$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$5!="0/1"&&$6!="0/1"&&$7!="0/1"' \
|awk '!($6=="1/1"&&$7=="1/1")'\
|awk '!($6=="0/0"&&$7=="0/0")'\
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> "$data"_merged_vcf_info_filtered_revised.txt;

done

## category_2 trio_vc_and_pedgree_phasing

for data in "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
HG002_phased_vcf=`echo "${!data}"`
mkdir $wkdir/$category_2/$data
cd $wkdir/$category_2/$data
# merge 3 vcf files 
bcftools merge $HG002_phased_vcf $HG003_phased_vcf $HG004_phased_vcf > "$data"_merged.vcf 

# extract  info from merged vcf 
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG002"\t"GQ_HG003"\t"GQ_HG004"\t"DP_HG002"\t"DP_HG003"\t"DP_HG004"\t"PS_HG002"\t"PS_HG003"\t"PS_HG004 > "$data"_merged_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' "$data"_merged.vcf >> "$data"_merged_vcf_info.txt

# note: no ./. called for HG002; however ./. will be intoduced for HG002 in the merged vcf file

# 1. remove ./. in HG002&HG003&HG004 & retain only biallelic snp 
# 2. remove unphased loci in HG003&HG004
# 3. 4. remove region where HG003==HG004==hom
# 5. 6. replace 0/0 1/1 with 0|0 1|1 
# 7.8.9.  replace 0/1 ./. in HG002 with .|.

head -1 "$data"_merged_vcf_info.txt>"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|awk '$5!="./."&&$6!="./."&&$7!="./."&&length($3)==1&&length($4)==1' \
|awk '$6!="0/1"&&$7!="0/1"' \
|awk '!($6=="1/1"&&$7=="1/1")'\
|awk '!($6=="0/0"&&$7=="0/0")'\
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'\
|sed 's#0/1#.|.#g'\
|sed 's#./.#.|.#g'>> "$data"_merged_vcf_info_filtered_revised.txt;

done