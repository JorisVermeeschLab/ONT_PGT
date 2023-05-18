#!/bin/bash

# chmod 755 4_get_summarized_file.sh
# nohup 4_get_summarized_file.sh>4_get_summarized_file.sh.log 2>&1 &
# r23i27n22
# 2900


wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_02_check_top_PGT_M_genes"

# gene info file path contain coorediates in GTF file
gene_info=$wkdir/gene_list_olga_with_info.txt

cd $wkdir

vc_phasing_dirs=`ls -d */|sed 's#/##'|awk 'BEGIN{ORS=" "}{print $1}'`

for vc_phasing_dir in $vc_phasing_dirs
do
all_data=`ls $vc_phasing_dir|awk 'BEGIN{ORS=" "}{print $1}'`

for data in $all_data
do 
path=$wkdir/$vc_phasing_dir/$data

# num_data_covered_gene=`cat $path/b_"$data"_trio_success.txt| grep -v CHROM|awk '{print $4}'|sort|uniq|wc -l`

# get stats 
echo -e "CHROM" "start" "end" "gene" "len_gene" "len_bp" "perc_gene" "len_snv" "max_matchperc_Ha" "max_matchperc_Hb" "class" "match_a" "match_b" "label"> $path/c_"$data"_sumarized_info.txt
awk 'NR==FNR{a[$4]=$3-$2+1;next}($4 in a){print $1,$2,$3,$4,a[$4],$5,$5/a[$4],$6,$22,$23,$(NF-2),$(NF-1),$NF,$1":"$2"_"$3}' $gene_info $path/b_"$data"_trio_success.txt|sort -k 14,14  >> $path/c_"$data"_sumarized_info.txt
done

done














