wkdir="/Users/miffy/OneDrive - KU Leuven/1_working directory/02_nanopore/03_human_PGT_families/ONT2/CNV/ONT2_E04/CNV"
sample=""
name=""
wd="$wkdir/$sample/CNV"
file=`ls $wkdir|grep HSL`

Rscript $wkdir/2_plot_cnv.R "$wkdir" "$sample" "$file" "$name"
