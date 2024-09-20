library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

data<-read.table(file=paste("a_gene_",sample,"_haplotype_match_stats.txt",sep=""),header = T)
# ensure dataframe format
data<-data.frame(data)

# ------------------A. autosomes and Pseudoautosomal regions ---------------------

data_sub_label<-data %>% filter(CHROM!="X")%>%
  rowwise()%>%
  mutate(match_a = paste(c("3_Ha","3_Hb")[c_across(perc_a_m_3_Ha:perc_a_m_3_Hb) == max_matchperc_Ha],collapse=",")) %>%
  mutate(match_b = paste(c("4_Ha","4_Hb")[c_across(perc_b_m_4_Ha:perc_b_m_4_Hb) == max_matchperc_Hb],collapse=",")) 

# ------------------B. chrx specific regions ---------------------

data_sub_uniq_x_label<-data %>% filter(CHROM=="X")%>%
              rowwise()%>%
              mutate(match_a = ".")%>%
              mutate(match_b = paste(c("4_Ha","4_Hb")[c_across(perc_b_m_4_Ha:perc_b_m_4_Hb) == max_matchperc_Hb_4_for_x],collapse=","))

# ------------------C. merge and output ---------------------
 
data_label<-rbind(data_sub_label,data_sub_uniq_x_label)

write.table(data_label,file=paste("b_",sample,"_label.txt",sep=""),row.names = F,quote = F)

