#!/bin/bash

# chmod 755 3_ped_check_gene_exon_cover.sh
# nohup 3_ped_check_gene_exon_cover.sh > 3_ped_check_gene_exon_cover.sh.log 2>&1 &
# 
# 
module load BEDTools

uniq_exon_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/GRCh38/nochr_noYM_uniq_exon_protein_coding.sorted.merged.bed"
uniq_gene_protein_coding_bed="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/GRCh38/nochr_noYM_uniq_gene_protein_coding.sorted.merged.bed"

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/GRCh38/1_check_GW_PGT_M/trio_vc_and_pedgree_phasing/check_protein_genes"

echo -e data"\t"len_snv_lim"\t"max_matchperc_lim"\t"exon_coverage"\t"gene_coverage>$wkdir/ped_gene_exon_coverage.txt

for data in "ped_mc_MDA_23x_LSK110" "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
cd $wkdir/$data

for len_snv_lim in 2 4 6
do
for max_matchperc_lim in 1 0.9 0.8 0.7 0.6 0
do

## 1.1 get phased block bed file --> GTF start positions -1 
cat b_"$data"_trio_success.txt|grep -v CHROM|sort -k1,1 -k2,2n|awk '$8>='''$len_snv_lim'''&& $24>='''$max_matchperc_lim'''&& $25>='''$max_matchperc_lim''' '|awk '{print $1"\t"$2-1"\t"$3}'|sort -k 1,1 -k2,2n|bedtools merge>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed

## 1.2 get stats for covered genes
echo -e CHROM"\t"start"\t"end"\t"gene_name"\t"gene_id"\t"gene_len_bp"\t"len_bp"\t"len_snv>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".coverd_info
cat b_"$data"_trio_success.txt|grep -v CHROM|sort -k1,1 -k2,2n|awk '$8>='''$len_snv_lim'''&& $24>='''$max_matchperc_lim'''&& $25>='''$max_matchperc_lim''' '|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}'>>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".coverd_info


## 2. get overlapped regions between uniq_exon_protein_coding_bed (GENCODE) and the sample phased block bed file

bedtools intersect -a c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed \
                             -b $uniq_exon_protein_coding_bed>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect.bed

## 3. get overlapped regions between uniq_gene_protein_coding_bed (GENCODE) and the sample phased block bed file      

bedtools intersect -a c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".bed \
                             -b $uniq_gene_protein_coding_bed>c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed

#awk '{sum_bp+=$3-$2}END{print sum_bp}' $uniq_exon_protein_coding_bed
# 108291662
sum_exon_protein_coding=108291662
#awk '{sum_bp+=$3-$2}END{print sum_bp}' $uniq_gene_protein_coding_bed
# 1301418355
sum_gene_protein_coding=1301418355

## 4. check gen/exon coverage 

exon_coverage_bed=c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".exon.intersect.bed
exon_covered_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $exon_coverage_bed`

gene_coverage_bed=c_"$data"_trio_success_lensnv_"$len_snv_lim"_max_matchperc_lim_"$max_matchperc_lim".gene.intersect.bed
gene_covered_success=`awk '{sum_bp+=$3-$2}END{print sum_bp}' $gene_coverage_bed`

exon_coverage=`echo -e $exon_covered_success"\t"$sum_exon_protein_coding|awk 'BEGIN { OFMT = "%.4f"} {print $1/$2}'`
gene_coverage=`echo -e $gene_covered_success"\t"$sum_gene_protein_coding|awk 'BEGIN { OFMT = "%.4f"} {print $1/$2}'`


echo -e $data"\t"$len_snv_lim"\t"$max_matchperc_lim"\t"$exon_coverage"\t"$gene_coverage>>$wkdir/ped_gene_exon_coverage.txt

rm c_*.bed
done
done
done



# --------------get max_matchperc stats----------------------
echo -e data"\t"min_max_matchperc_Ha"\t"min_max_matchperc_Hb"\t"perc_blocks_max_matchperc_Ha_Hb_1>$wkdir/ped_min_max_matchperc.txt
for data in "ped_mc_MDA_23x_LSK110" "ped_sc_MDA_24x_LSK110" "ped_mc_MDA_24x_LSK114" "ped_sc_MDA_22x_LSK114" "ped_bulk_no_24x_LSK114"
do
cd $wkdir/$data

# Force a type conversion, using addition of zero
min_max_matchperc_Ha=`cat b_"$data"_trio_success.txt|grep -v CHROM|awk 'BEGIN{a=1}{if ($24<0+a) a=$24}END{print a}'`
min_max_matchperc_Hb=`cat b_"$data"_trio_success.txt|grep -v CHROM|awk 'BEGIN{a=1}{if ($25<0+a) a=$25}END{print a}'`

nrow=`cat b_"$data"_trio_success.txt|grep -v CHROM|wc -l`

perc_blocks_max_matchperc_Ha_Hb_1=`cat b_"$data"_trio_success.txt|grep -v CHROM|awk '$24==1&&$25==1'|awk -v nrow=$nrow 'END{print NR/nrow}'`


echo -e $data"\t"$min_max_matchperc_Ha"\t"$min_max_matchperc_Hb"\t"$perc_blocks_max_matchperc_Ha_Hb_1 >>$wkdir/ped_min_max_matchperc.txt

done 