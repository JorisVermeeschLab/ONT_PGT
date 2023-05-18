#!/bin/bash

# chmod 755 4_gene_exon_covered_for_bed_regions_data_overlap_with_ref.sh
# nohup 4_gene_exon_covered_for_bed_regions_data_overlap_with_ref.sh>4_gene_exon_covered_for_bed_regions_data_overlap_with_ref.sh.log 2>&1 &
# r23i27n24
# 26840

module load BEDTools

uniq_exon_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GENCODE/gtf_bed_files/nochr_uniq_exon_protein_coding.sorted.merged.bed"
uniq_gene_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GENCODE/gtf_bed_files/nochr_uniq_gene_protein_coding.sorted.merged.bed"

GIAB_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_01_resource_GIAB_HG002_ref/4.2.1/1_HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.sorted.bed"

ref="GIAB_4_2_1_hifiasm_v11_phasetransfer"

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_1_2_test_explore_PGT_M_merge_whole_vcf"

echo -e data"\t"len_snv_lim"\t"max_matchperc_lim"\t"exon_coverage"\t"gene_coverage"\t"perc_overlap_vs_ref_exon"\t"perc_overlap_vs_ref_gene"\t"perc_overlap_vs_data_exon"\t"perc_overlap_vs_data_gene>$wkdir/gene_exon_coverage_compare_with_ref_GIAB_bed_regions.txt

## 1. intersect the lensnv_2_max_matchperc_lim_0 ref bed files with GIAB bed region
cd $wkdir/$ref

bedtools intersect -a c_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.exon.intersect.bed -b $GIAB_bed > d_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.exon.intersect_GIAB_bed_region.bed
bedtools intersect -a c_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.gene.intersect.bed -b $GIAB_bed > d_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.gene.intersect_GIAB_bed_region.bed

compare_exon="$wkdir/$ref/d_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.exon.intersect_GIAB_bed_region.bed"
compare_gene="$wkdir/$ref/d_"$ref"_trio_success_lensnv_2_max_matchperc_lim_0.gene.intersect_GIAB_bed_region.bed"

for data in "bulk_no_10x_LSK110" "sc_PTA_10x_LSK110" "mc_PTA_4x_LSK110" "sc_MDA_24x_LSK110" "mc_MDA_23x_LSK110" "sc_MDA_22x_LSK114" "mc_MDA_24x_LSK114" "bulk_no_24x_LSK114"
do
cd $wkdir/$data

for len_snv_lim in 2 4 6
do
for max_matchperc_lim in 1 0.9 0.8 0.7 0.6 0
do
# 2.intersect data covered gene/exon regions with GIAB bed regions 

# exon
bedtools intersect -a c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect.bed \
                             -b $GIAB_bed > d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect_GIAB_bed_region.bed

# gene
bedtools intersect -a c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed \
                             -b $GIAB_bed > d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect_GIAB_bed_region.bed

## 4.get overlap for data vs ref 

## get overlapped regions exon
bedtools intersect -a d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect_GIAB_bed_region.bed \
                             -b $compare_exon > e_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect_GIAB_bed_region_overlap_data_ref.bed
## get overlapped regions gene  
bedtools intersect -a d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect_GIAB_bed_region.bed \
                             -b $compare_gene > e_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect_GIAB_bed_region_overlap_data_ref.bed


# ref covered gene/exon regions whitin the GIAB bed region (lensnv_2_max_matchperc_lim_0)
sum_bed_ref_exon_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $compare_exon`
sum_bed_fef_gene_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $compare_gene`

# data covered gene/exon regions whitin the GIAB bed region 
sum_bed_data_exon_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect_GIAB_bed_region.bed`
sum_bed_data_gene_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' d_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect_GIAB_bed_region.bed`

# ref & data both covered gene/exon regions whitin the GIAB bed region 
sum_bed_overlap_success_exon=`awk '{sum_bp+=$3-$2}END{print sum_bp}' e_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect_GIAB_bed_region_overlap_data_ref.bed`
sum_bed_overlap_success_gene=`awk '{sum_bp+=$3-$2}END{print sum_bp}' e_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect_GIAB_bed_region_overlap_data_ref.bed`



perc_overlap_vs_ref_exon=`echo -e $sum_bed_overlap_success_exon"\t"$sum_bed_ref_exon_success|awk 'BEGIN {OFMT = "%.4f"}{print $1/$2}'`
perc_overlap_vs_ref_gene=`echo -e $sum_bed_overlap_success_gene"\t"$sum_bed_fef_gene_success|awk 'BEGIN {OFMT = "%.4f"}{print $1/$2}'`

perc_overlap_vs_data_exon=`echo -e $sum_bed_overlap_success_exon"\t"$sum_bed_data_exon_success|awk 'BEGIN {OFMT = "%.4f"}{print $1/$2}'`
perc_overlap_vs_data_gene=`echo -e $sum_bed_overlap_success_gene"\t"$sum_bed_data_gene_success|awk 'BEGIN {OFMT = "%.4f"}{print $1/$2}'`

## 5. check gen/exon coverage data 

#awk '{sum_bp+=$3-$2}END{print sum_bp}' $uniq_exon_protein_coding_bed
# 108522233
sum_exon_protein_coding=108522233
#awk '{sum_bp+=$3-$2}END{print sum_bp}' $uniq_gene_protein_coding_bed
# 1304978080
sum_gene_protein_coding=1304978080

exon_coverage_bed=c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect.bed
exon_covered_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $exon_coverage_bed`

gene_coverage_bed=c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed
gene_covered_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $gene_coverage_bed`

exon_coverage=`echo -e $exon_covered_success"\t"$sum_exon_protein_coding|awk 'BEGIN { OFMT = "%.4f"} {print $1/$2}'`
gene_coverage=`echo -e $gene_covered_success"\t"$sum_gene_protein_coding|awk 'BEGIN { OFMT = "%.4f"} {print $1/$2}'`

echo -e $data"\t"$len_snv_lim"\t"$max_matchperc_lim"\t"$exon_coverage"\t"$gene_coverage"\t"$perc_overlap_vs_ref_exon"\t"$perc_overlap_vs_ref_gene"\t"$perc_overlap_vs_data_exon"\t"$perc_overlap_vs_data_gene>>$wkdir/gene_exon_coverage_compare_with_ref_GIAB_bed_regions.txt


done
done
done