library(dplyr)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

to_check<-read.table(paste0("d_",sample,"_to_check_hap_infer.txt"),header=T)%>%
            mutate(informativity=case_when(GT_HG003_called=="0/0"&GT_HG004_called=="0/0"~"noninformative",
                                           GT_HG003_called=="1/1"&GT_HG004_called=="1/1"~"noninformative",
                                           GT_HG003_called=="0/1"|GT_HG004_called=="0/1"~"noninformative",
                                           GT_HG002_called=="noninformative"~"noninformative",
                                           TRUE~"informative"))


to_check$CHROM<-as.character(to_check$CHROM)

# find ID for each loci in to_check
# ID_loci<-read.table(paste0("a_",sample,"_gt_sep_ps_imputed.txt"),header=T)%>%

file="/lustre1/project/stg_00019/research/yan/nanopore_data/01_01_HG002_sc_mc_bulk/03_01_PGT_M/4_hg38_use_patmat_match_FMr10/2_E_ped_phasing/0_get_unfiltered_block_info/a_FM_gt_sep_ps_imputed.txt"
ID_loci<-read.table(file,header=T)%>%
select(CHROM,POS,REF,ALT,ID)

to_check_ID<-left_join(to_check,ID_loci,join_by(CHROM,POS,REF,ALT))

# # get inferred GT
hap_label<-read.table(paste0("b_",sample,"_label.txt"),header=T)

hp_PGT_performance<-left_join(to_check_ID, hap_label, join_by(CHROM,POS>=start,POS<=end,ID))%>% 
	     mutate(inferred_GT=case_when(informativity=="informative"&(match_a=="3_Ha"|match_a=="3_Ha,3_Hb")&(match_b=="4_Ha"|match_b=="4_Ha,4_Hb")~paste0(substr(GT_HG003_called,1,1),"|",substr(GT_HG004_called,1,1)),
	                                        informativity=="informative"&match_a=="3_Hb"&match_b=="4_Hb"~paste0(substr(GT_HG003_called,3,3),"|",substr(GT_HG004_called,3,3)),
	                                        informativity=="informative"&(match_a=="3_Ha"|match_a=="3_Ha,3_Hb")&(match_b=="4_Hb")~paste0(substr(GT_HG003_called,1,1),"|",substr(GT_HG004_called,3,3)),
	                                        informativity=="informative"&(match_a=="3_Hb")&(match_b=="4_Ha"|match_b=="4_Ha,4_Hb")~paste0(substr(GT_HG003_called,3,3),"|",substr(GT_HG004_called,1,1)),
                                          TRUE~"no"))%>%
  select(CHROM,POS,GT_HG002_ref,inferred_GT,informativity,pattern,type)

write.table(hp_PGT_performance,paste0("d_",sample,"_haplotype_compare.txt"),quote=F,col.names=T,row.names=F)


