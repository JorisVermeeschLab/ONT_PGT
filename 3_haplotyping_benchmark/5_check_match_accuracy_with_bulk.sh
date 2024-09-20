#!/bin/bash

# chmod 755 5_check_match_accuracy_with_bulk.sh
# ./5_check_match_accuracy_with_bulk.sh

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing"
bulk_label="$wkdir/ped_bulk/b_ped_bulk_label.txt"
n_block_bulk_label=`cat $bulk_label|grep -v CHROM|wc -l`

for data in "ped_mc" "ped_sc"
do 
cd $wkdir/$data
data_label=b_"$data"_label.txt
n_block_data_label=`cat $data_label|grep -v CHROM|wc -l`

file1=$bulk_label
file2=$data_label
#get colunns from both files for same regions   
echo -e ID"\t"match_a"\t"match_b"\t"match_a_bulk"\t"match_b_bulk> e_"$data"_compare_with_bulk.txt
awk 'NR==FNR{a[$4]=$(NF-1)"\t"$NF;next}($4 in a){print $4"\t"$(NF-1)"\t"$NF"\t"a[$4]}' $file1 $file2|grep -v ID >> e_"$data"_compare_with_bulk.txt

# get stats for data 
n_compared=`cat e_"$data"_compare_with_bulk.txt|grep -v ID|wc -l`
n_match_bulk=`cat e_"$data"_compare_with_bulk.txt|awk '$2==$4&&$3==$5'|wc -l`
perc_matched_of_compared=`echo -e $n_match_bulk"\t"$n_compared|awk 'BEGIN { OFMT = "%.4f"}{print $1/$2}'`

echo -e n_block_bulk_label"\t"$n_block_bulk_label"\n"\
n_block_data_label"\t"$n_block_data_label"\n"\
n_compared"\t"$n_compared"\n"\
n_match_bulk"\t"$n_match_bulk"\n"\
perc_matched_of_compared"\t"$perc_matched_of_compared>e_accuracy_"$data"_stats.txt;

rm e_"$data"_compare_with_bulk.txt

done



