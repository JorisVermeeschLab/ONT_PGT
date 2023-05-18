
cd /lustre1/project/stg_00019/research/yan/nanopore_data/03_01_2_check_top_PGT_M_genes

# changed some names in the list to correct version
genelist="gene_list_olga.txt"
encode_list="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GENCODE/gtf_bed_files/gene_protein_coding.txt"

awk 'NR==FNR{a[$6]=$1"\t"$2"\t"$3"\t"$6;next}($0 in a){print a[$0]}' $encode_list $genelist>gene_list_olga_with_info.txt

# remove chr in gene_list_olga_with_info.txt

awk 'NR==FNR{a[$6]=$1"\t"$2"\t"$3"\t"$6;next}!($0 in a){print $0}' $encode_list $genelist>gene_list_olga_no_info.txt