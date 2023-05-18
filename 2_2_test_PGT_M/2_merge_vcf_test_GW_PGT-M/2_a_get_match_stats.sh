# nohup 2_a_get_match_stats.sh > 2_a_get_match_stats.sh.log 2>&1 &
# r23i27n24
# 21108

module load R/4.0.2-foss-2018a-bare


wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_1_2_test_explore_PGT_M_merge_whole_vcf"
cd $wkdir

for data in "bulk_no_10x_LSK110" "sc_PTA_10x_LSK110" "mc_PTA_4x_LSK110" "sc_MDA_24x_LSK110" "mc_MDA_23x_LSK110" "sc_MDA_22x_LSK114" "mc_MDA_24x_LSK114" "bulk_no_24x_LSK114"
do
Rscript a_get_match_stats.R "$wkdir/$data" $data
done

