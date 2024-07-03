#!/bin/bash

# chmod 755 4_check_match_accuracy.sh
# ./4_check_match_accuracy.sh

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/GRCh38/1_check_GW_PGT_M/trio_vc_and_pedgree_phasing/check_protein_genes"
bulk_trio=$wkdir/ped_bulk_no_24x_LSK114/b_ped_bulk_no_24x_LSK114_trio_success.txt

# compare sc/mc dataset with ref (bulk trio) 

for data in "ped_mc_MDA_23x_LSK110" "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114"
do 
path=$wkdir/$data
file=$path/b_"$data"_trio_success.txt

#get columns from both files for same regions   
echo -e CHROM"\t"start"\t"end"\t"gene_name"\t"gene_id"\t"class_bulk"\t"match_a_bulk"\t"match_b_bulk"\t"class"\t"match_a"\t"match_b> $path/d_"$data"_compare_info.txt
awk 'NR==FNR{a[$1$2$3$4$5]=$(NF-2)"\t"$(NF-1)"\t"$NF;next}($1$2$3$4$5 in a){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$(NF-2)"\t"$(NF-1),$NF,a[$1$2$3$4$5]}' $bulk_trio $file|grep -v CHROM >> $path/d_"$data"_compare_info.txt

# get stats
n_block_ref=`cat $bulk_trio|grep -v CHROM|wc -l`
n_block_data=`cat $file|grep -v CHROM|wc -l`
n_block_merged=`cat $path/d_"$data"_compare_info.txt|grep -v CHROM|wc -l`
n_block_merged_matched=`cat $path/d_"$data"_compare_info.txt|grep -v CHROM|awk '($(NF)==$(NF-3))&&($(NF-1)==$(NF-4))&&($(NF-2)==$(NF-5))'|wc -l`
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
>$path/d_accuracy_"$data"_stats.txt;

done



