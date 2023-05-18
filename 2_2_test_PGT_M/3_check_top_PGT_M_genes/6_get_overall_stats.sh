#!/bin/bash

# chmod 755 
# nohup .log 2>&1 &
# r23i27n24
# 

bed="/lustre1/project/stg_00019/research/yan/nanopore_data/03_02_check_top_PGT_M_genes/gene_list.bed"
bp_genes=`cat $bed|awk '{sum+=($3-$2);}END{print sum;}'`

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_02_check_top_PGT_M_genes"

cd $wkdir
echo -e vc_phasing"\t"data"\t"num_genes_covered"\t"perc_num_genes_covered"\t"perc_bp_genes_covered"\t"perc_merged_bulk_trio"\t"perc_merged_data"\t"perc_matched_vs_merged>top_genes_stat.txt

vc_phasing_dirs=`ls -d */|sed 's#/##'|awk 'BEGIN{ORS=" "}{print $1}'`

for vc_phasing_dir in $vc_phasing_dirs
do
all_data=`ls $vc_phasing_dir|awk 'BEGIN{ORS=" "}{print $1}'`

for data in $all_data
do 
path=$wkdir/$vc_phasing_dir/$data

file_summarize=$path/c_"$data"_sumarized_info.txt
file_stat=$path/d_merge_bulk_trio_"$data"_stats.txt

num_genes_covered=`cat $file_summarize|grep -v CHROM|awk '{print $4}'|sort|uniq|wc -l`
perc_num_genes_covered=`echo -e $num_genes_covered"\t"201|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`

bp_genes_covered=`cat $file_summarize|grep -v CHROM|awk '{sum+=$6;}END{print sum;}'`
perc_bp_genes_covered=`echo -e $bp_genes_covered"\t"$bp_genes|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`

perc_merged_bulk_trio=`cat $file_stat|grep perc_merged_bulk_trio|head -1|cut -f 2`
perc_merged_data=`cat $file_stat|grep perc_merged_"$data"|cut -f 2`
perc_matched_vs_merged=`cat $file_stat|grep perc_matched_vs_merged|cut -f 2`

echo -e $vc_phasing_dir"\t"$data"\t"$num_genes_covered"\t"$perc_num_genes_covered"\t"$perc_bp_genes_covered"\t"$perc_merged_bulk_trio"\t"$perc_merged_data"\t"$perc_matched_vs_merged>>top_genes_stat.txt
done
done