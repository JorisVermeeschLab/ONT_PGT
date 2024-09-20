#!/bin/bash
# chmod 755 6_block_stats.sh
# nohup 6_block_stats.sh > 6_block_stats.sh.log 2>&1 &
#  
# conda activate "/lustre1/project/stg_00019/research/yan/conda_env/R_4.3"

for data in ped_bulk ped_mc ped_sc
do
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing/$data"
Rscript 6_block_stats.R $wkdir $data;
done;