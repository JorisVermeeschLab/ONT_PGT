library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

# read in merged_vcf_info_filtered_revised_rm_small_block.txt
file<-paste(sample,"_merged_vcf_info_filtered_revised.txt",sep="")
data_raw<-read.table(file,header = T)

# first filter based on only parents
# for chrx specific region: HG004 retain only phased het loci  & HG003 change  GT & GQ 
# PAR1,PAR12: remove either unphased & both hom loci  
data_MIC_labeled<-data_raw %>% 
       mutate(CHROM_new = case_when(CHROM=="X"&(POS>=10001&POS<=2781479) ~ "PAR1", 
                                   CHROM=="X"&(POS>=155701383&POS<=156030895) ~ "PAR2", 
                                   TRUE ~ CHROM)) %>% 
       filter(!(CHROM_new=="X"&(GT_HG004=="0/1"|GT_HG004=="0|0"|GT_HG004=="1|1")))  %>% 
       mutate(GT_HG003=case_when(CHROM_new=="X" ~ "NC|NC",
                                    TRUE ~ GT_HG003)) %>%
       mutate(GQ_HG003=case_when(CHROM_new=="X" ~ "100",
                                    TRUE ~ GQ_HG003)) %>%
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_HG003=="0/1"|GT_HG004=="0/1"))) %>% 
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_HG003=="1|1"&GT_HG004=="1|1")))%>% 
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_HG003=="0|0"&GT_HG004=="0|0")))%>%
      mutate(type = case_when(
      CHROM_new!="X"&GT_HG003=="1|1" & (GT_HG004=="0|1"|GT_HG004=="1|0") & GT_HG002=="0|0" ~ "MIC",
      CHROM_new!="X"&GT_HG004=="1|1" & (GT_HG003=="0|1"|GT_HG003=="1|0") & GT_HG002=="0|0" ~ "MIC",
      CHROM_new!="X"&GT_HG003=="0|0" & (GT_HG004=="0|1"|GT_HG004=="1|0") & GT_HG002=="1|1" ~ "MIC",
      CHROM_new!="X"&GT_HG004=="0|0" & (GT_HG003=="0|1"|GT_HG003=="1|0") & GT_HG002=="1|1" ~ "MIC",
      CHROM_new!="X"&GT_HG003=="1|1" & GT_HG004=="0|0" & (GT_HG002=="1|1"|GT_HG002=="0|0") ~ "MIC_ADO",
      CHROM_new!="X"&GT_HG004=="1|1" & GT_HG003=="0|0" & (GT_HG002=="1|1"|GT_HG002=="0|0") ~ "MIC_ADO", 
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

# change PS on chrX specific regions for HG003 --> make the 3 regions distinct. 
# same paternal PS -- > blocks are separated by only maternal PS for chrx specific regions.

data<-data_MIC_labeled %>% 
      mutate(PS_HG003=case_when(CHROM_new=="X"&POS<10001 ~ "1",
                          CHROM_new=="X"&POS>2781479&POS<155701383 ~ "2781480",
                          CHROM_new=="X"&POS>156030895 ~ "156030896",
                          TRUE ~ PS_HG003)) 


 # replace . in PS_** columns with value from the surrounding line
 # filter out more son's loci and MIC loci
data_ps_imputed<-data %>%
  group_by(CHROM_new) %>%
  mutate(PS_HG003 = as.numeric(na_if(PS_HG003, '.'))) %>%
  mutate(PS_HG004 = as.numeric(na_if(PS_HG004, '.'))) %>%
  mutate(PS_HG002 = as.numeric(na_if(PS_HG002, '.'))) %>%
  fill(PS_HG003,PS_HG004,PS_HG002,.direction = "downup") %>% 
  ungroup() %>% 
  filter(type=="normal") %>%  
  select(-type) %>%   
  filter(GT_HG002!="0/1") %>%   
  filter(!(CHROM_new=="X"&(GT_HG002=="1|0"|GT_HG002=="0|1")))

# group by CHROM_new & PS; gave each group a unique ID; PS_HG002 is the same within a chr; PS_HG003+PS_HG004 enough to label a block
# numbering to ensure length<1MB 
# seperate GT 
data_gt_sep_ps_imputed<-data_ps_imputed %>%                   
  mutate(numbering = POS %/% 1000000+1) %>%
  mutate(ID=paste(CHROM_new,PS_HG003,PS_HG004,numbering,sep="_")) %>%                      
  separate(GT_HG003,c("HG003_Ha","HG003_Hb"),"\\|")  %>%
  separate(GT_HG004,c("HG004_Ha","HG004_Hb"),"\\|")  %>%
  separate(GT_HG002,c("HG002_Ha","HG002_Hb"),"\\|")

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
              seq_2_Ha=paste(HG002_Ha,collapse = ""),
              seq_2_Hb=paste(HG002_Hb,collapse = ""),
              seq_3_Ha=paste(HG003_Ha,collapse = ""),
              seq_3_Hb=paste(HG003_Hb,collapse = ""),
              seq_4_Ha=paste(HG004_Ha,collapse = ""),
              seq_4_Hb=paste(HG004_Hb,collapse = ""),
              a_m_3_Ha=sum(HG003_Ha==HG002_Ha), #  get the number of matched loci for each parental-son haplotype pair
              a_m_3_Hb=sum(HG003_Hb==HG002_Ha),
              b_m_4_Ha=sum(HG004_Ha==HG002_Hb),
              b_m_4_Hb=sum(HG004_Hb==HG002_Hb)) %>%
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


