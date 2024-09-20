#!/bin/bash

# chmod 755 2_a_get_match_stats.sh 
# nohup 2_a_get_match_stats.sh > 2_a_get_match_stats.sh.log 2>&1 &
# r23i27n22  1052453
# conda activate "/lustre1/project/stg_00019/research/yan/conda_env/R_4.3"

for data in ONT2_E04 ONT2_E06 ONT2_E20
do
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/PGT_M_snv_patmat_match_new_model/$data"
Rscript 2_a_get_match_stats.R $wkdir $data;
done