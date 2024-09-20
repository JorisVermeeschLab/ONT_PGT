#!/bin/bash
# chmod 755 select_informative_loci.sh
# nohup select_informative_loci.sh>select_informative_loci.sh.log 2>&1 &
#  

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BCFtools/1.9-foss-2018a

SAMPLE=""

# path to PGT-M analysis dir
wkdir=""
cd $wkdir

# default varinat calling results of the parents and the offspring
F_filtered_vcf="xx/6_Clair3/filtered.merge_output.vcf.gz"
M_filtered_vcf="xx/6_Clair3/filtered.merge_output.vcf.gz"
E_vcf="xx/6_Clair3/merge_output.vcf.gz"  

zcat $E_vcf|grep '^#'> $SAMPLE.filtered.merge_output.vcf
zcat $E_vcf|grep -v '#'|awk -F: '$(NF-2)>2'>> $SAMPLE.filtered.merge_output.vcf
bgzip $SAMPLE.filtered.merge_output.vcf

E_filtered_vcf=$SAMPLE.filtered.merge_output.vcf.gz
tabix $E_filtered_vcf

# merge vcf files 
bcftools merge --missing-to-ref $E_filtered_vcf $F_filtered_vcf $M_filtered_vcf|grep -v RefCall> "$SAMPLE"_merged.vcf 

# get vcf info, including AF(frequency for an alternate allele)
# retain loci where the parental showed different homozygous genotypes
# 00 11 or 11 00
# extract  info from merged vcf 
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_E"\t"GT_F"\t"GT_M"\t"AF_E"\t"AF_F"\t"AF_M"\t"DP_E"\t"DP_F"\t"DP_M > diff_ref_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%AF][\t%DP]\n' "$SAMPLE"_merged.vcf \
|awk 'length($3)==1&&length($4)==1'\
|awk '$1!="Y"' \
|awk '($6=="0/0"&&$7=="1/1")||($6=="1/1"&&$7=="0/0")' \
|sed 's#0/0#0#g'\
|sed 's#1/1#1#g'\
|sed 's#\t\.#\t0#g'>> diff_ref_vcf_info.txt

# extract info from merged vcf (default) for aneuoloid chr 

echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_E"\t"GT_F"\t"GT_M"\t"AF_E"\t"AF_F"\t"AF_M"\t"DP_E"\t"DP_F"\t"DP_M > check_origion.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%AF][\t%DP]\n' "$SAMPLE"_merged.vcf \
|awk 'length($3)==1&&length($4)==1'\
|awk '($6=="0/1"&&($7=="1/1"||$7=="0/0"))||($7=="0/1"&&($6=="1/1"||$6=="0/0"))'>> check_origion.txt