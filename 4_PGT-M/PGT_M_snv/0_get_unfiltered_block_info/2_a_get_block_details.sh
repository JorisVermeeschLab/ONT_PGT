#!/bin/bash

# chmod 755 2_a_get_block_details.sh 
# nohup 2_a_get_block_details.sh > 2_a_get_block_details.sh.log 2>&1 &
# i28l05  
# conda activate "/lustre1/project/stg_00019/research/yan/conda_env/R_4.3"


wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/PGT_M_snv_patmat_match_new_model/0_get_unfiltered_block_info"
Rscript 2_a_get_block_details.R $wkdir; 
