library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

# read in merged_vcf_info_filtered_revised_rm_small_block.txt
file<-paste(sample,"_merged_vcf_info_filtered_revised.txt",sep="")
data_raw<-read.table(file,header = T)

# first filter based on only parents
# for chrx specific region: M retain only phased het loci  & F change  GT & GQ 
# PAR1,PAR12: remove either unphased & both hom loci  
data_MIC_labeled<-data_raw %>% 
       mutate(CHROM_new = case_when(CHROM=="X"&(POS>=10001&POS<=2781479) ~ "PAR1", 
                                   CHROM=="X"&(POS>=155701383&POS<=156030895) ~ "PAR2", 
                                   TRUE ~ CHROM)) %>% 
       filter(!(CHROM_new=="X"&(GT_M=="0/1"|GT_M=="0|0"|GT_M=="1|1")))  %>% 
       mutate(GT_F=case_when(CHROM_new=="X" ~ "NC|NC",
                                    TRUE ~ GT_F)) %>%
       mutate(GQ_F=case_when(CHROM_new=="X" ~ "100",
                                    TRUE ~ GQ_F)) %>%
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_F=="0/1"|GT_M=="0/1"))) %>% 
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_F=="1|1"&GT_M=="1|1")))%>% 
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_F=="0|0"&GT_M=="0|0")))%>%
      mutate(type = case_when(
      CHROM_new!="X"&GT_F=="1|1" & (GT_M=="0|1"|GT_M=="1|0") & GT_E=="0|0" ~ "MIC",
      CHROM_new!="X"&GT_M=="1|1" & (GT_F=="0|1"|GT_F=="1|0") & GT_E=="0|0" ~ "MIC",
      CHROM_new!="X"&GT_F=="0|0" & (GT_M=="0|1"|GT_M=="1|0") & GT_E=="1|1" ~ "MIC",
      CHROM_new!="X"&GT_M=="0|0" & (GT_F=="0|1"|GT_F=="1|0") & GT_E=="1|1" ~ "MIC",
      CHROM_new!="X"&GT_F=="1|1" & GT_M=="0|0" & (GT_E=="1|1"|GT_E=="0|0") ~ "MIC_ADO",
      CHROM_new!="X"&GT_M=="1|1" & GT_F=="0|0" & (GT_E=="1|1"|GT_E=="0|0") ~ "MIC_ADO", 
                          TRUE ~ "normal")) 
MIC<-data_MIC_labeled%>%
        filter(type!="normal")

# stats_MIC for autosomes and PAR1 PAR2
stats_MIC<-data_MIC_labeled%>%
            filter(CHROM_new!="X")%>%
            group_by(CHROM_new)%>%
            summarize(
		  	    num_loci=n(),
		  	    num_MIC=sum(type!="normal"),       
		  	    num_ADO=sum(type=="MIC_ADO"))

# overall stats_MIC for autosomes and PAR1 PAR2
overall_stats_MIC<-stats_MIC %>% summarize(
				   CHROM_new="autosomes_PARs",
				   num_loci=sum(num_loci),
				   num_MIC=sum(num_MIC),
				   num_ADO=sum(num_ADO))

final_stats_MIC<-rbind(stats_MIC,overall_stats_MIC)

final_stats_MIC_with_perc<-final_stats_MIC %>% 
mutate(
perc_MIC=num_MIC/num_loci,
perc_ADO=num_ADO/num_loci)

# output stats_MIC and MIC
write.table(final_stats_MIC_with_perc,file=paste("a_",sample,"_stats_MIC_autosomes_PARs.txt",sep=""),row.names = F,quote = F)
write.table(MIC,file=paste("a_",sample,"_MIC_autosomes_PARs.txt",sep=""),row.names = F,quote = F)

# change PS on chrX specific regions for F --> make the 3 regions distinct. 
# same paternal PS -- > blocks are separated by only maternal PS for chrx specific regions.

data<-data_MIC_labeled %>% 
      mutate(PS_F=case_when(CHROM_new=="X"&POS<10001 ~ "1",
                          CHROM_new=="X"&POS>2781479&POS<155701383 ~ "2781480",
                          CHROM_new=="X"&POS>156030895 ~ "156030896",
                          TRUE ~ PS_F)) 


 # replace . in PS_** columns with value from the surrounding line
 # filter out more son's loci and MIC loci
data_ps_imputed<-data %>%
  group_by(CHROM_new) %>%
  mutate(PS_F = as.numeric(na_if(PS_F, '.'))) %>%
  mutate(PS_M = as.numeric(na_if(PS_M, '.'))) %>%
  mutate(PS_E = as.numeric(na_if(PS_E, '.'))) %>%
  fill(PS_F,PS_M,PS_E,.direction = "downup") %>% 
  ungroup() %>% 
  filter(type=="normal") %>%  
  select(-type) %>%   
  filter(GT_E!="0/1") %>%   
  filter(!(CHROM_new=="X"&(GT_E=="1|0"|GT_E=="0|1")))

# group by CHROM_new & PS; gave each group a unique ID; PS_E is the same within a chr; PS_F+PS_M enough to label a block
# numbering to ensure length<1MB 
# seperate GT 
data_gt_sep_ps_imputed<-data_ps_imputed %>%                   
  mutate(numbering = POS %/% 1000000+1) %>%
  mutate(ID=paste(CHROM_new,PS_F,PS_M,numbering,sep="_")) %>%                      
  separate(GT_F,c("F_Ha","F_Hb"),"\\|")  %>%
  separate(GT_M,c("M_Ha","M_Hb"),"\\|")  %>%
  separate(GT_E,c("E_Ha","E_Hb"),"\\|")

write.table(data_gt_sep_ps_imputed,file=paste("a_",sample,"_gt_sep_ps_imputed.txt",sep=""),row.names = F,quote = F)

# 1. group by ID 
# 2. get stats for each group 
# retain blocks with >=2 SNVs
# start and end positions are from VCF file (1 based)

  trio_PS_stats<-
    data_gt_sep_ps_imputed %>%
    group_by(ID) %>%
    summarise(CHROM=unique(CHROM_new),
              start=min(POS),
              end=max(POS),
              len_bp=max(POS)-min(POS)+1,
              len_snv= n(), # get the number of row per group
              seq_2_Ha=paste(E_Ha,collapse = ""),
              seq_2_Hb=paste(E_Hb,collapse = ""),
              seq_3_Ha=paste(F_Ha,collapse = ""),
              seq_3_Hb=paste(F_Hb,collapse = ""),
              seq_4_Ha=paste(M_Ha,collapse = ""),
              seq_4_Hb=paste(M_Hb,collapse = ""),
              a_m_3_Ha=sum(F_Ha==E_Ha), #  get the number of matched loci for each parental-son haplotype pair
              a_m_3_Hb=sum(F_Hb==E_Ha),
              b_m_4_Ha=sum(M_Ha==E_Hb),
              b_m_4_Hb=sum(M_Hb==E_Hb)) %>%
              filter(len_snv>=2) %>% 
              arrange(CHROM,start) %>% 
              ungroup()

clean_trio_PS_stats<-trio_PS_stats %>%
  select(CHROM,start,end,ID,len_bp,len_snv,a_m_3_Ha,a_m_3_Hb,b_m_4_Ha,b_m_4_Hb)

# consider unique regions on chrx 

final_trio_match_stats<-clean_trio_PS_stats %>%
  rowwise() %>%
  mutate(max_matchperc_Ha=max(across(contains("a_m_")))/len_snv, max_matchperc_Hb=max(across(contains("b_m_")))/len_snv) %>%
  mutate(perc_a_m_3_Ha=a_m_3_Ha/len_snv) %>% 
  mutate(perc_a_m_3_Hb=a_m_3_Hb/len_snv) %>% 
  mutate(perc_b_m_4_Ha=b_m_4_Ha/len_snv) %>% 
  mutate(perc_b_m_4_Hb=b_m_4_Hb/len_snv)%>%
  mutate(e_match_Ha_3=sum(across(contains("perc_a_m_3_")) == max_matchperc_Ha)) %>%
  mutate(e_match_Hb_4=sum(across(contains("perc_b_m_4_")) == max_matchperc_Hb))  %>%
  mutate(max_matchperc_Hb_4_for_x=max(across(contains("b_m_4")))/len_snv) %>%
  mutate(e_match_Hb_4_for_x=sum(across(contains("perc_b_m_4")) == max_matchperc_Hb_4_for_x))

write.table(final_trio_match_stats,file=paste("a_gene_",sample,"_haplotype_match_stats.txt",sep=""),row.names = F,quote = F)


