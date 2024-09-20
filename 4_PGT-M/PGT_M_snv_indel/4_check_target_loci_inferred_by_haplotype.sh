#!/bin/bash
# chmod 755 4_check_target_loci_inferred_by_haplotype.sh
# nohup 4_check_target_loci_inferred_by_haplotype.sh > 4_check_target_loci_inferred_by_haplotype.sh.log 2>&1 &
# 
# conda activate "/lustre1/project/stg_00019/research/yan/conda_env/R_4.3"

for sample in ONT2_E04 ONT2_E06 ONT2_E20
do
wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/PGT_M_snv_indel_patmat_match_new_model/$sample"
cd $wkdir
# var_locus_1="1 209618002"  # del
# var_locus_2="1 209628898"  # SNV

# get unit ID
ID_1=`cat a_"$sample"_gt_sep_ps_imputed.txt|grep ^"1 209618002"|awk '{print $NF}'`
ID_2=`cat a_"$sample"_gt_sep_ps_imputed.txt|grep ^"1 209628898"|awk '{print $NF}'`

# out put all iSNP in the unit
head -1 a_"$sample"_gt_sep_ps_imputed.txt|awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$NF}'> d_variant_unit_detail_1.txt 
cat a_"$sample"_gt_sep_ps_imputed.txt|awk -v var=$ID_1 '$NF==var{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$NF}' >> d_variant_unit_detail_1.txt

head -1 a_"$sample"_gt_sep_ps_imputed.txt|awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$NF}'> d_variant_unit_detail_2.txt 
cat a_"$sample"_gt_sep_ps_imputed.txt|awk -v var=$ID_2 '$NF==var{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$NF}' >> d_variant_unit_detail_2.txt

# get match stats
echo -e CHROM start end len_bp len_snv H1_inferred H2_inferred H1_m_F1 H1_m_F2 H2_m_M1 H2_m_M2 ID>d_variant_unit_match_stats_1.txt 
cat b_"$sample"_label.txt|grep $ID_1|awk '{print $1,$2,$3,$5,$6,$(NF-1),$NF,$7,$8,$9,$10,$4}'>>d_variant_unit_match_stats_1.txt

echo -e CHROM start end len_bp len_snv H1_inferred H2_inferred H1_m_F1 H1_m_F2 H2_m_M1 H2_m_M2 ID>d_variant_unit_match_stats_2.txt 
cat b_"$sample"_label.txt|grep $ID_2|awk '{print $1,$2,$3,$5,$6,$(NF-1),$NF,$7,$8,$9,$10,$4}'>>d_variant_unit_match_stats_2.txt

# out put all SNP in the unit
info="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/PGT_M_snv_indel_patmat_match_new_model/0_get_unfiltered_block_info/a_FM_gt_sep_ps_imputed.txt"
head -1 $info|awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$NF}'> d_variant_unit_detail_all_1.txt 
cat $info|awk -v var=$ID_1 '$NF==var{print $1,$2,$3,$4,$5,$6,$7,$8,$NF}' >> d_variant_unit_detail_all_1.txt

head -1 $info|awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$NF}'> d_variant_unit_detail_all_2.txt 
cat $info|awk -v var=$ID_2 '$NF==var{print $1,$2,$3,$4,$5,$6,$7,$8,$NF}' >> d_variant_unit_detail_all_2.txt

done