#!/bin/bash

# chmod 755 2_a_get_match_stats.sh 
# nohup 2_a_get_match_stats.sh > 2_a_get_match_stats.sh.log 2>&1 &
# 
# 

module load R/4.0.2-foss-2018a-bare

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/2_check_top_PGT_M_genes"

#independent_vc_and_independent_phasing
wkdir_ind="$wkdir/independent_vc_and_independent_phasing"

for data in "independent_sc_PTA_10x_LSK110" "independent_mc_PTA_4x_LSK110" "independent_sc_MDA_24x_LSK110" "independent_mc_MDA_23x_LSK110" "independent_mc_MDA_24x_LSK114" "independent_sc_MDA_22x_LSK114" "independent_bulk_no_24x_LSK114" 
do
Rscript a_get_match_stats.R "$wkdir_ind/$data" $data;
done;

# trio_vc_and_pedgree_phasing
wkdir_ped="$wkdir/trio_vc_and_pedgree_phasing"

for data in "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_23x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
Rscript a_get_match_stats.R "$wkdir_ped/$data" $data;
done
