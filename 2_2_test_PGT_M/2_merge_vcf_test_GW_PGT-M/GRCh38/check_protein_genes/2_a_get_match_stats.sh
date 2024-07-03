#!/bin/bash

# chmod 755 2_a_get_match_stats.sh 
# nohup 2_a_get_match_stats.sh > 2_a_get_match_stats.sh.log 2>&1 &
# 31320

module load R/4.0.2-foss-2018a-bare

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/GRCh38/1_check_GW_PGT_M/trio_vc_and_pedgree_phasing/check_protein_genes"

for data in "ped_mc_MDA_23x_LSK110" "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
Rscript 2_a_get_match_stats.R "$wkdir/$data" $data;
done;
