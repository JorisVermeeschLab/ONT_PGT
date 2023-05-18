library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

file<-paste(sample,"_merged_vcf_info_filtered_revised.txt",sep="")
data_raw<-read.table(file,header = T)

data<-data_raw %>%
filter(GT_HG002!=".|.")

 # replace . in PS_** columns with value from the surrounding line
data_ps_imputed<-data %>%
  group_by(CHROM) %>%
  mutate(PS_HG003 = as.numeric(na_if(PS_HG003, '.'))) %>%
  mutate(PS_HG004 = as.numeric(na_if(PS_HG004, '.'))) %>%
  mutate(PS_HG002 = as.numeric(na_if(PS_HG002, '.'))) %>%
  fill(PS_HG003,PS_HG004,PS_HG002,.direction = "downup") # #default direction down

# merge PS from the trio to a single colunm  CHROM is also added 
data_ps_imputed$PS_merged<-paste(data_ps_imputed$CHROM,data_ps_imputed$PS_HG003,data_ps_imputed$PS_HG004,data_ps_imputed$PS_HG002)

# seperate GT in the tri to sepqrate allele
data_gt_sep_ps_imputed<-data_ps_imputed %>%
  separate(GT_HG003,c("HG003_Ha","HG003_Hb"),"\\|")  %>%
  separate(GT_HG004,c("HG004_Ha","HG004_Hb"),"\\|")  %>%
  separate(GT_HG002,c("HG002_Ha","HG002_Hb"),"\\|")

# group data_ps_imputed by PS_merged
# get stats for each PS_merged group
  trio_PS_stats<-
    data_gt_sep_ps_imputed %>%
    group_by(PS_merged) %>%
    summarise(CHROM=unique(CHROM),
              start=min(POS),
              end=max(POS),
              len_bp=max(POS)-min(POS)+1,
              len_snv= n(), # nrow by group
              seq_2_Ha=paste(HG002_Ha,collapse = ""),
              seq_2_Hb=paste(HG002_Hb,collapse = ""),
              seq_3_Ha=paste(HG003_Ha,collapse = ""),
              seq_3_Hb=paste(HG003_Hb,collapse = ""),
              seq_4_Ha=paste(HG004_Ha,collapse = ""),
              seq_4_Hb=paste(HG004_Hb,collapse = ""),
              a_m_3_Ha=sum(HG003_Ha==HG002_Ha),
              a_m_3_Hb=sum(HG003_Hb==HG002_Ha),
              a_m_4_Ha=sum(HG004_Ha==HG002_Ha),
              a_m_4_Hb=sum(HG004_Hb==HG002_Ha),
              b_m_3_Ha=sum(HG003_Ha==HG002_Hb),
              b_m_3_Hb=sum(HG003_Hb==HG002_Hb),
              b_m_4_Ha=sum(HG004_Ha==HG002_Hb),
              b_m_4_Hb=sum(HG004_Hb==HG002_Hb))  %>%
              filter(len_snv>=2) %>% 
              arrange(CHROM,start)
           

trio_PS_stats<-trio_PS_stats %>%  ungroup()

# check 1. num unique hp for parental hp & 2. repeat times for each parental hp

trio_PS_hp_stats<-trio_PS_stats %>% 
  rowwise() %>%
   mutate(rep_3_Ha = sum(c_across(seq_3_Ha:seq_4_Hb) == seq_3_Ha),
          rep_3_Hb = sum(c_across(seq_3_Ha:seq_4_Hb) == seq_3_Hb),
          rep_4_Ha = sum(c_across(seq_3_Ha:seq_4_Hb) == seq_4_Ha),
          rep_4_Hb = sum(c_across(seq_3_Ha:seq_4_Hb) == seq_4_Hb),
          unique_hp_FP = n_distinct(c_across(seq_3_Ha:seq_4_Hb)),
          unique_hp_3 = n_distinct(c_across(seq_3_Ha:seq_3_Hb)),
          unique_hp_4 = n_distinct(c_across(seq_4_Ha:seq_4_Hb))) %>%
          ungroup()

clean_trio_PS_stats<-trio_PS_hp_stats %>%
  select(CHROM,start,end,len_bp,len_snv,a_m_3_Ha,a_m_3_Hb,a_m_4_Ha,a_m_4_Hb,b_m_3_Ha,b_m_3_Hb,b_m_4_Ha,b_m_4_Hb,unique_hp_FP,unique_hp_3,unique_hp_4,rep_3_Ha,rep_3_Hb,rep_4_Ha,rep_4_Hb)

final_trio_match_stats<-clean_trio_PS_stats %>%
  rowwise() %>%
  mutate(max_matchperc_Ha=max(across(contains("a_m_")))/len_snv, max_matchperc_Hb=max(across(contains("b_m_")))/len_snv) %>%
  mutate(perc_a_m_3_Ha=a_m_3_Ha/len_snv) %>% 
  mutate(perc_a_m_3_Hb=a_m_3_Hb/len_snv) %>% 
  mutate(perc_a_m_4_Ha=a_m_4_Ha/len_snv) %>% 
  mutate(perc_a_m_4_Hb=a_m_4_Hb/len_snv) %>% 
  mutate(perc_b_m_3_Ha=b_m_3_Ha/len_snv) %>% 
  mutate(perc_b_m_3_Hb=b_m_3_Hb/len_snv) %>% 
  mutate(perc_b_m_4_Ha=b_m_4_Ha/len_snv) %>% 
  mutate(perc_b_m_4_Hb=b_m_4_Hb/len_snv)%>%
  mutate(e_match_Ha=sum(across(contains("perc_a_m_")) == max_matchperc_Ha)) %>%
  mutate(e_match_Hb=sum(across(contains("perc_b_m_")) == max_matchperc_Hb)) %>%
  mutate(e_match_Ha_3=sum(across(contains("perc_a_m_3_")) == max_matchperc_Ha)) %>%
  mutate(e_match_Ha_4=sum(across(contains("perc_a_m_4_")) == max_matchperc_Ha)) %>%
  mutate(e_match_Hb_3=sum(across(contains("perc_b_m_3_")) == max_matchperc_Hb)) %>%
  mutate(e_match_Hb_4=sum(across(contains("perc_b_m_4_")) == max_matchperc_Hb)) 

write.table(final_trio_match_stats,file=paste("a_",sample,"_haplotype_match_stats.txt",sep=""),row.names = F,quote = F)


