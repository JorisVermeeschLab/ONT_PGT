#!/bin/bash
# chmod 755 1_merge_vcf_and_extract_required_info.sh
# nohup 1_merge_vcf_and_extract_required_info.sh>1_merge_vcf_and_extract_required_info.sh.log 2>&1 &
# r23i27n22 1018406

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BCFtools/1.9-foss-2018a

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/PGT_M_snv_patmat_match_new_model"

F_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/ONT2_F/6_Clair3/whatshap_phasing_GQ2/phased.filtered.merge_output.vcf.gz"
M_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/ONT2_M/6_Clair3/whatshap_phasing_GQ2/phased.filtered.merge_output.vcf.gz"

for data in ONT2_E04 ONT2_E06 ONT2_E20
do

E_ped_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/ONT2/$data/7_Clair3_trio_newmodel/whatshap_phasing/$data.phased.vcf.gz"

outdir=$wkdir/$data
mkdir -p $outdir
cd $outdir

# merge vcf files 
bcftools merge --missing-to-ref $E_ped_phased_vcf $F_phased_vcf $M_phased_vcf > "$data"_merged.vcf 

# extract info from merged vcf; remove loci with only 0/0
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_E"\t"GT_F"\t"GT_M"\t"GQ_E"\t"GQ_F"\t"GQ_M"\t"DP_E"\t"DP_F"\t"DP_M"\t"PS_E"\t"PS_F"\t"PS_M > "$data"_merged_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' "$data"_merged.vcf \
|awk '!($5=="0/0"&&$6=="0/0"&&$7=="0/0")'>> "$data"_merged_vcf_info.txt

# filter the parental loci 

head -1 "$data"_merged_vcf_info.txt> $outdir/"$data"_merged_vcf_info_filtered_revised.txt

cat "$data"_merged_vcf_info.txt\
|grep -v "CHROM" \
|awk 'length($3)==1&&length($4)==1' \
|awk '$1!="Y"' \
|awk '!($1!="X"&&($6=="0/1"||$7=="0/1"))' \
|awk '!($1!="X"&&($6=="1/1"&&$7=="1/1"))' \
|awk '!($1!="X"&&($6=="0/0"&&$7=="0/0"))' \
|sed 's#0/0#0|0#g'\
|sed 's#1/1#1|1#g'>> $outdir/"$data"_merged_vcf_info_filtered_revised.txt;

rm "$data"_merged.vcf 
done
