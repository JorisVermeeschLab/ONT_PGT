#!/bin/bash
# chmod 755 1_get_single_merged_FMvcf_file.sh
# nohup 1_get_single_merged_FMvcf_file.sh>1_get_single_merged_FMvcf_file.sh.log 2>&1 &
#  

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BCFtools/1.9-foss-2018a

F_phased_vcf="xx/phased.filtered.merge_output.vcf.gz"
M_phased_vcf="xx/phased.filtered.merge_output.vcf.gz"


wkdir=""

# merge FM phased vcf files 
bcftools merge --missing-to-ref $F_phased_vcf $M_phased_vcf > parents.vcf 

# extract info from merged vcf 
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG003"\t"GQ_HG004"\t"PS_HG003"\t"PS_HG004 > parents_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%PS]\n' parents.vcf  >> parents_vcf_info.txt

# filter the parental loci 

head -1 parents_vcf_info.txt> parents_vcf_info_filtered_revised.txt

cat parents_vcf_info.txt\
|grep -v "CHROM" \
|awk 'length($3)==1&&length($4)==1' \
|awk '$1!="Y"' \
|awk '!($1!="X"&&($5=="0/1"||$6=="0/1"))' \
|awk '!($1!="X"&&($5=="1/1"&&$6=="1/1"))' \
|awk '!($1!="X"&&($5=="0/0"&&$6=="0/0"))' \
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> parents_vcf_info_filtered_revised.txt;

rm parents.vcf 
