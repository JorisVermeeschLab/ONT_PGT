
cd /lustre1/project/stg_00019/research/yan/nanopore_data/03_01_PGT_M/2_GRCh38_CHM13v1.0_with_refcall_with_GQ_filtering/CHM13v1.0/2_check_top_PGT_M_genes

# GFF3 gene info dir
info_dir="/lustre1/project/stg_00019/research/yan/nanopore_data/00_02_resource_GTF/T2T-CHM13v1.0"

# changed some names in the list to correct version
genelist="gene_list_olga.txt"
anno_list="$info_dir/chm13.draft_v1.0_protein_gene_info.txt"

awk 'NR==FNR{a[$6]=$1"\t"$3"\t"$4"\t"$6;next}($0 in a){print a[$0]}' $anno_list $genelist>gene_list_olga_with_info.txt

awk 'NR==FNR{a[$6]=$1"\t"$3"\t"$4"\t"$6;next}!($0 in a){print $0}' $anno_list $genelist>genes_not_found.txt