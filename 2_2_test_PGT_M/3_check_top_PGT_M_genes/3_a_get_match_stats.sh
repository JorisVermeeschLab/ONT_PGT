# chmod 755 3_a_get_match_stats.sh
# nohup 3_a_get_match_stats.sh > 3_a_get_match_stats.sh.log 2>&1 &
# r23i27n22
# 35009

module load R/4.0.2-foss-2018a-bare

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_02_check_top_PGT_M_genes"

for dir in "independent_vc" "trio_variant_calling" "trio_vc_phsing" "trio_vc_phsing_merge_ONT_trio"
do
data_all=`ls $wkdir/$dir|awk 'BEGIN{ORS=" "}{print}'`
for data in $data_all
do
R_wkdir=$wkdir/$dir/$data
Rscript a_get_match_stats.R "$R_wkdir" $data;
done
done


