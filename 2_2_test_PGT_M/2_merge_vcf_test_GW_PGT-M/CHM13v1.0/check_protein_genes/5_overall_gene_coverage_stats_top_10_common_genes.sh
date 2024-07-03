#!/bin/bash

# chmod 755 5_overall_gene_coverage_stats_top_10_common_genes.sh
# nohup 5_overall_gene_coverage_stats_top_10_common_genes.sh 5_overall_gene_coverage_stats_top_10_common_genes.sh.log 2>&1 &

# /lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/T2T-CHM13v1.0/gene_list_with_info.txt
# cat gene_list_with_info.txt|grep -v seq|wc -l 
# 19956 unique gene ID on chr1-22-x 
# num_gene_CHM13v1.0=19956

# /lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/GRCh38/gene_list_with_info.txt
# cat gene_list_with_info.txt|grep -v seq|wc -l 
# 19954 unique gene ID on chr1-22-x
# num_gene_GRCh38=19954

module load R/4.0.2-foss-2018a-bare

num_gene_ID_protein_coding=19956

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/1_check_GW_PGT_M/trio_vc_and_pedgree_phasing/check_protein_genes"

# for data in "ped_mc_MDA_23x_LSK110" "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
for data in "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
Rscript 5_overall_gene_coverage_stats_top_10_common_genes.R "$wkdir/$data" $data $num_gene_ID_protein_coding;
done