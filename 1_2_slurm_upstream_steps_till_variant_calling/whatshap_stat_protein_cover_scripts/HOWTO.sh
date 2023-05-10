#!/bin/bash

# chmod 755 HOWTO.sh 1_whatshap_stats_protein_cover.sh 2_check_covered_exon_gene_stat.R
# nohup HOWTO.sh > whatshap_stat_protein_cover.log 2>&1 &

wkdir="xx_variant_calling"
SAMPLE="xx" 
wdR="$wkdir/protein_gene_exon_cover"
./1_whatshap_stats_protein_cover.sh $wkdir $SAMPLE
Rscript 2_check_covered_exon_gene_stat.R $SAMPLE $wdR