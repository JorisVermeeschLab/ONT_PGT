# #!/bin/bash

# chmod 755 2_a_get_match_stats.sh a_get_match_stats.R
# nohup 2_a_get_match_stats.sh > 2_a_get_match_stats.sh.log 2>&1 &
# r23i27n22
# 2217

module load R/4.0.2-foss-2018a-bare

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/1_check_GW_PGT_M"
cd $wkdir
#independent_vc_and_independent_phasing
wkdir_ind="$wkdir/independent_vc_and_independent_phasing"

for data in "independent_mc_MDA_24x_LSK114" "independent_sc_MDA_22x_LSK114" "independent_bulk_no_24x_LSK114" 
do
Rscript 2_a_get_match_stats.R "$wkdir_ind/$data" $data;
done

# trio_vc_and_pedgree_phasing
wkdir_ped="$wkdir/trio_vc_and_pedgree_phasing"

for data in "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
Rscript 2_a_get_match_stats.R "$wkdir_ped/$data" $data;
done
