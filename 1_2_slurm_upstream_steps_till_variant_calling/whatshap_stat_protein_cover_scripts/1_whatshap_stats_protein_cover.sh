#!/bin/bash -l
#PBS -l walltime=6:00:00
#PBS -l mem=10gb
#PBS -l nodes=1:ppn=2
#PBS -M yan.zhao1@student.kuleuven.be
#PBS -m eab
#PBS -A lp_joris_vermeesch_c1


uniq_exon_protein_coding_bed=/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GENCODE/gtf_bed_files/nochr_uniq_exon_protein_coding.bed
uniq_gene_protein_coding_bed=/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GENCODE/gtf_bed_files/nochr_uniq_gene_protein_coding.bed

# SAMPLE and wkdir to be changed 

phased_vcf=$1/phased_merge_output.vcf.gz
out_dir=$1/protein_gene_exon_cover
mkdir -p $out_dir
cd $out_dir

#check vcfstats
/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/rtg-tools-3.12.1/rtg vcfstats $phased_vcf>$2.vcfstats.txt

# create a GTF file from a phased VCF file that describes the haplotype blocks.
# get phasing_statistics
whatshap stats --gtf=$2.phased_merge_output.gtf $phased_vcf>$2.phasing_statistics.txt

# cat $uniq_exon_protein_coding_bed |wc -l
# 502253
# cat $uniq_exon_protein_coding_bed|awk '$2!=$3'|wc -l
# 502253
# cat $uniq_gene_protein_coding_bed |wc -l
# 20028
# cat $uniq_gene_protein_coding_bed|awk '$2!=$3'|wc -l
# 20028

# change the temp folder by using the TMPDIR environment variable
export TMPDIR=$VSC_SCRATCH

## 1.get phased block bed file from gtf file  --> GTF start positions -1 
cat $2.phased_merge_output.gtf|awk '{print $1"\t"$4-1"\t"$5}'>$2.phased_block.bed

## 2. get overlapped regions between uniq_exon_protein_coding_bed (GENCODE) and the sample phased block bed file
module load BEDTools

bedtools intersect -a $2.phased_block.bed \
                             -b $uniq_exon_protein_coding_bed>$2.exon.intersect.bed

## 3. get overlapped regions between uniq_geneprotein_coding_bed (GENCODE) and the sample phased block bed file      
bedtools intersect -a $2.phased_block.bed \
                             -b $uniq_gene_protein_coding_bed>$2.gene.intersect.bed               

## 4. deduplication of the overlapped regions
cat $2.exon.intersect.bed|sort|uniq>$2.exon.intersect.uniq.bed
cat $2.gene.intersect.bed|sort|uniq>$2.gene.intersect.uniq.bed

## 5. get the number of overlapped regions that cover whole genes 
#  -d, --repeated  only print duplicate lines, one for each group
cat $uniq_exon_protein_coding_bed $2.exon.intersect.uniq.bed|sort|uniq -d>$2.exon.cover.bed
n_exon=`wc -l $2.exon.cover.bed`
echo "complete protein coding exons covered" $n_exon "all protein coding exons 502253"

cat $uniq_gene_protein_coding_bed $2.gene.intersect.uniq.bed|sort|uniq -d>$2.gene.cover.bed
n_gene=`wc -l $2.gene.cover.bed`
echo "complete protein coding genes covered" $n_gene "all protein coding genes 20028"

rm $2.phased_block.bed $2.exon.intersect.bed $2.gene.intersect.bed $2.exon.intersect.uniq.bed $2.gene.intersect.uniq.bed