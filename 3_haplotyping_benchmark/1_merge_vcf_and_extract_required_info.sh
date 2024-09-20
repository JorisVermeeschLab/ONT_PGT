#!/bin/bash
# chmod 755 1_merge_vcf_and_extract_required_info.sh
# nohup 1_merge_vcf_and_extract_required_info.sh>1_merge_vcf_and_extract_required_info.sh.log 2>&1 &
# i28l05 2824575

module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load BCFtools/1.9-foss-2018a

wkdir="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing"

HG003_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/giab_202305/HG003_F/6_Clair3/whatshap_phasing_GQ2/phased.filtered.merge_output.vcf.gz"
HG004_phased_vcf="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/giab_202305/HG004_M/6_Clair3/whatshap_phasing_GQ2/phased.filtered.merge_output.vcf.gz"

ped_mc="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/01_02_data_LSK114_HG002/1_mc_MDA_24x_LSK114/7_variant_calling_trio_parents_newkit/whatshap_phasing/mc_MDA_24x_LSK114.phased.vcf.gz"
ped_sc="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/01_02_data_LSK114_HG002/2_sc_MDA_22x_LSK114/7_variant_calling_trio_parents_newkit/whatshap_phasing/sc_MDA_22x_LSK114.phased.vcf.gz"
ped_bulk="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/01_02_data_LSK114_HG002/3_bulk_no_24x_LSK114/7_variant_calling_trio_parents_newkit/whatshap_phasing/bulk_no_24x_LSK114.phased.vcf.gz"

for data in "ped_mc" "ped_sc" "ped_bulk"

do

outdir=$wkdir/$data
mkdir -p $outdir
cd $outdir

HG002_phased_vcf=`echo "${!data}"`

# merge vcf files 
bcftools merge --missing-to-ref $HG002_phased_vcf $HG003_phased_vcf $HG004_phased_vcf > "$data"_merged.vcf 

# extract info from merged vcf; remove loci with only 0/0
echo -e CHROM"\t"POS"\t"REF"\t"ALT"\t"GT_HG002"\t"GT_HG003"\t"GT_HG004"\t"GQ_HG002"\t"GQ_HG003"\t"GQ_HG004"\t"DP_HG002"\t"DP_HG003"\t"DP_HG004"\t"PS_HG002"\t"PS_HG003"\t"PS_HG004 > "$data"_merged_vcf_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%DP][\t%PS]\n' "$data"_merged.vcf \
|awk '!($5=="0/0"&&$6=="0/0"&&$7=="0/0")'>> "$data"_merged_vcf_info.txt

# get stats for informative loci on autosomes
cat "$data"_merged_vcf_info.txt \
|grep -v "CHROM" \
|awk 'length($3)==1&&length($4)==1' \
|awk '$1!="Y"&&$1!="X"' \
|awk '!($6=="0/1"||$7=="0/1")' \
|awk '!($6=="1/1"&&$7=="1/1")' \
|awk '!($6=="0/0"&&$7=="0/0")' \
|awk '{print $5,$6,$7}' \
|sort|uniq -c>"$data"_phase_GT_stats.txt

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