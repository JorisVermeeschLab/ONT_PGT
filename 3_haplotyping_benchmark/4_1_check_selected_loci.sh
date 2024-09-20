#!/bin/bash

# chmod 755 4_1_check_selected_loci.sh
# nohup 4_1_check_selected_loci.sh>4_1_check_selected_loci.sh.log 2>&1 &
#   

loci_dir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/0_select_snp_indel_GIAB"
selected_info="$loci_dir/PGT_M_candidates_only_snv.txt"

ind_vcf_dir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/1_E_ind_phasing"
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing"

for data in ped_bulk ped_mc ped_sc
do

ind_data=`echo $data|sed 's/ped/ind/'`
cd $wkdir/$data

# -------------evaluate accuracy of original called genotype -------------
total=`cat $selected_info|grep -v CHROM|wc -l`

file1=$ind_vcf_dir/$ind_data/"$ind_data"_merged_vcf_info.txt 
file2=$selected_info

# get comparison info
echo -e CHROM"\t"POS"\t"GT_HG002_ref"\t"GT_HG002_called"\t"pattern"\t"type>d_"$data"_direct_mutation_comparison.txt
# a. selected ref called in data
awk 'NR==FNR{a[_$1_$2_$3_$4_]=$5; next} ((_$1_$2_$3_$4_) in a) {print $1"\t"$2"\t"$5"\t"a[_$1_$2_$3_$4_]"\t"$8"\t"$9}' $file1 $file2 \
|grep -v CHROM>>d_"$data"_direct_mutation_comparison.txt
# b. selected ref nocall in data-->impute to 0/0
awk 'NR==FNR{a[_$1_$2_$3_$4_]=$5; next} !((_$1_$2_$3_$4_) in a) {print $1"\t"$2"\t"$5"\t""0/0""\t"$8"\t"$9}' $file1 $file2 \
|grep -v CHROM>>d_"$data"_direct_mutation_comparison.txt

# -------------get d_loci_to_check_hap_infer;-------------

file1="$data"_merged_vcf_info.txt 
file2=$selected_info

# get comparison info
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002_ref"\t"GT_HG002_called"\t"GT_HG003_called"\t"GT_HG004_called"\t"pattern"\t"type>d_"$data"_to_check_hap_infer.txt
# a. selected ref available in data
awk 'NR==FNR{a[_$1_$2_$3_$4_]=$5" "$6" "$7; next} ((_$1_$2_$3_$4_) in a) {print $1"\t"$2"\t"$3"\t"$4"\t"$5,a[_$1_$2_$3_$4_],$8,$9}' $file1 $file2\
|grep -v CHROM>>d_"$data"_to_check_hap_infer.txt
# b. selected ref noninformative(unavailable) in data
awk 'NR==FNR{a[_$1_$2_$3_$4_]=$5" "$6" "$7; next} !((_$1_$2_$3_$4_) in a) {print $1"\t"$2"\t"$3"\t"$4"\t"$5,"noninformative","noninformative","noninformative",$8,$9}' $file1 $file2\
|grep -v CHROM>>d_"$data"_to_check_hap_infer.txt
done