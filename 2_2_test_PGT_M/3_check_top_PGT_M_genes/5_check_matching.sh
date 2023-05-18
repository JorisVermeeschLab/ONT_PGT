#!/bin/bash

# chmod 755 5_check_matching.sh
# nohup 5_check_matching.sh>5_check_matching.sh.log 2>&1 &
# r23i27n24
# 4625


wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_02_check_top_PGT_M_genes"
ref=$wkdir/independent_vc/GIAB_4_2_1_hifiasm_v11_phasetransfer/c_GIAB_4_2_1_hifiasm_v11_phasetransfer_sumarized_info.txt
bulk_trio=$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/c_bulk_no_24x_LSK114_sumarized_info.txt

# # check similarity ref vs bulk_trio

# #merge 
# echo -e "CHROM" "start" "end" "gene" "len_gene" "len_bp" "perc_gene" "len_snv_ref" "len_snv_data" "max_matchperc_Ha_ref" "max_matchperc_Hb_ref" "max_matchperc_Ha_data" "max_matchperc_Hb_data" "class_ref" "class_data" "match_a_ref" "match_b_ref" "match_a_data" "match_b_data">$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/merge_ref_bulk_trio.txt
# join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $ref $bulk_trio|grep -v CHROM>>$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/merge_ref_bulk_trio.txt
# # get unmatched item in merged file
# echo -e "CHROM" "start" "end" "gene" "len_gene" "len_bp" "perc_gene" "len_snv_ref" "len_snv_data" "max_matchperc_Ha_ref" "max_matchperc_Hb_ref" "max_matchperc_Ha_data" "max_matchperc_Hb_data" "class_ref" "class_data" "match_a_ref" "match_b_ref" "match_a_data" "match_b_data">$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/merge_ref_bulk_trio_unmatched.txt
# join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $ref $bulk_trio|grep -v CHROM|awk '!($(NF-1)==$(NF-2)&&$(NF)==$(NF-3)||$(NF-1)==$(NF-3)&&$(NF)==$(NF-2))'>>$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/merge_ref_bulk_trio_unmatched.txt

# # get stats
# n_block_ref=`cat $ref|wc -l`
# n_block_data=`cat $bulk_trio|wc -l`
# n_block_merged=`join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $ref $bulk_trio|wc -l`
# n_block_merged_matched=`join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $ref $bulk_trio|awk '$(NF-1)==$(NF-2)&&$(NF)==$(NF-3)||$(NF-1)==$(NF-3)&&$(NF)==$(NF-2)'|wc -l`
# perc_merged_ref=`echo -e $n_block_merged"\t"$n_block_ref|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`
# perc_merged_data=`echo -e $n_block_merged"\t"$n_block_data|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`
# perc_matched_vs_merged=`echo -e $n_block_merged_matched"\t"$n_block_merged|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`

# echo -e n_block_ref"\t"$n_block_ref"\n"\
# n_block_bulk_trio"\t"$n_block_bulk_trio"\n"\
# n_block_merged"\t"$n_block_merged"\n"\
# n_block_merged_matched"\t"$n_block_merged_matched"\n"\
# perc_merged_ref"\t"$perc_merged_ref"\n"\
# perc_merged_bulk_trio"\t"$perc_merged_data"\n"\
# perc_matched_vs_merged"\t"$perc_matched_vs_merged"\n"\
# >$wkdir/trio_vc_phsing/bulk_no_24x_LSK114/merge_ref_bulk_trio_stats.txt

# compare each dataset with ref and bulk trio 
cd $wkdir

vc_phasing_dirs=`ls -d */|sed 's#/##'|awk 'BEGIN{ORS=" "}{print $1}'`

for vc_phasing_dir in $vc_phasing_dirs
do
all_data=`ls $vc_phasing_dir|awk 'BEGIN{ORS=" "}{print $1}'`

for data in $all_data
do 
path=$wkdir/$vc_phasing_dir/$data

file=$path/c_"$data"_sumarized_info.txt

# use bulk_trio as ref; merge with bulk_trio 

#merge 
echo -e "CHROM" "start" "end" "gene" "len_gene" "len_bp" "perc_gene" "len_snv_ref" "len_snv_data" "max_matchperc_Ha_ref" "max_matchperc_Hb_ref" "max_matchperc_Ha_data" "max_matchperc_Hb_data" "class_ref" "class_data" "match_a_ref" "match_b_ref" "match_a_data" "match_b_data">$path/d_merge_bulk_trio_"$data".txt
join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $bulk_trio $file|grep -v CHROM >> $path/d_merge_bulk_trio_"$data".txt

# get unmatched item in merged file
echo -e "CHROM" "start" "end" "gene" "len_gene" "len_bp" "perc_gene" "len_snv_ref" "len_snv_data" "max_matchperc_Ha_ref" "max_matchperc_Hb_ref" "max_matchperc_Ha_data" "max_matchperc_Hb_data" "class_ref" "class_data" "match_a_ref" "match_b_ref" "match_a_data" "match_b_data">$path/d_merge_bulk_trio_"$data"_unmatched.txt
join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $bulk_trio $file|grep -v CHROM|awk '!($(NF-1)==$(NF-2)&&$(NF)==$(NF-3)||$(NF-1)==$(NF-3)&&$(NF)==$(NF-2))'>>$path/d_merge_bulk_trio_"$data"_unmatched.txt

# get stats
n_block_ref=`cat $bulk_trio|wc -l`
n_block_data=`cat $file|wc -l`
n_block_merged=`join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $bulk_trio $file|wc -l`
n_block_merged_matched=`join -j 14 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,1.9,1.10,2.9,2.10,1.11,2.11,1.12,1.13,2.12,2.13 $bulk_trio $file|awk '$(NF-1)==$(NF-2)&&$(NF)==$(NF-3)||$(NF-1)==$(NF-3)&&$(NF)==$(NF-2)'|wc -l`
perc_merged_ref=`echo -e $n_block_merged"\t"$n_block_ref|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`
perc_merged_data=`echo -e $n_block_merged"\t"$n_block_data|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`
perc_matched_vs_merged=`echo -e $n_block_merged_matched"\t"$n_block_merged|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`

echo -e n_block_bulk_trio"\t"$n_block_ref"\n"\
n_block_$data"\t"$n_block_data"\n"\
n_block_merged"\t"$n_block_merged"\n"\
n_block_merged_matched"\t"$n_block_merged_matched"\n"\
perc_merged_bulk_trio"\t"$perc_merged_ref"\n"\
perc_merged_$data"\t"$perc_merged_data"\n"\
perc_matched_vs_merged"\t"$perc_matched_vs_merged"\n"\
>$path/d_merge_bulk_trio_"$data"_stats.txt;

done
done