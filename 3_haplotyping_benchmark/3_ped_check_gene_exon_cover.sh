#!/bin/bash

# chmod 755 3_ped_check_gene_exon_cover.sh
# nohup 3_ped_check_gene_exon_cover.sh > 3_ped_check_gene_exon_cover.sh.log 2>&1 &
# i28l05 
# 

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BEDTools

uniq_exon_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/GRCh38/nochr_noYM_uniq_exon_protein_coding.sorted.merged.bed"
uniq_gene_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/GRCh38/nochr_noYM_uniq_gene_protein_coding.sorted.merged.bed"


for data in ped_bulk ped_mc ped_sc
do
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing/$data"
cd $wkdir
echo -e data"\t"len_snv_lim"\t"max_matchperc_lim"\t"gene_coverage>c_gene_exon_coverage.txt

for len_snv_lim in 2 4 6 
do

for max_matchperc_lim in 1 0.9 0.8 0.7 0.6 0.5
do
cat b_"$data"_label.txt|grep -v CHROM|awk '$1!="X"&&$6>='''$len_snv_lim'''&& $11>='''$max_matchperc_lim'''&& $12>='''$max_matchperc_lim''' '|awk '{print $1"\t"$2-1"\t"$3}'>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed
cat b_"$data"_label.txt|grep -v CHROM|awk '$1=="X"&&$6>='''$len_snv_lim'''&& $19>='''$max_matchperc_lim''' '|awk '{print $1"\t"$2-1"\t"$3}'>>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed


## get overlapped regions between uniq_gene_protein_coding_bed (GENCODE) and the sample phased block bed file      

bedtools intersect -a c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed \
                             -b $uniq_gene_protein_coding_bed>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed

## 4. check gene coverage 

#awk '{sum_bp+=$3-$2}END{print sum_bp}' $uniq_gene_protein_coding_bed
# 1301418355
sum_gene_protein_coding=1301418355

gene_covered=`awk '{sum_bp+=$3-$2}END{print sum_bp}' c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed`

gene_coverage=`echo -e $gene_covered"\t"$sum_gene_protein_coding|awk 'BEGIN { OFMT = "%.4f"} {print $1/$2}'`

echo -e $data"\t"$len_snv_lim"\t"$max_matchperc_lim"\t"$gene_coverage>>c_gene_exon_coverage.txt;

rm *.bed 

done;
done;
done;