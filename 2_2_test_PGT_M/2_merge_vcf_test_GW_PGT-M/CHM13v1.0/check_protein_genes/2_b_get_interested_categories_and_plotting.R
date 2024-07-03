library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

data_for_selection<-read.table(file=paste("a_gene_",sample,"_haplotype_match_stats.txt",sep=""),header = T)

# ensure dataframe format
data_for_selection<-data.frame(data_for_selection)

#-----; do not know real unique_hp for the son. infer from the match pattern with parental hp--------
# consider rep_3_Ha rep_3_Hb rep_4_Ha rep_4_Hb 

## 1. unique_hp_FP 2

# FM2_003_2_004_2

#FM2_003_2_004_2_002_1  Ha Hb match to same parental hp
FM2_003_2_004_2_002_1<-data_for_selection %>%
  filter(unique_hp_FP == 2,unique_hp_3==2,unique_hp_4==2,e_match_Ha==2,e_match_Hb==2,(perc_a_m_3_Ha==max_matchperc_Ha&perc_b_m_3_Ha==max_matchperc_Hb)|(perc_a_m_3_Hb==max_matchperc_Ha&perc_b_m_3_Hb==max_matchperc_Hb))
#FM2_003_2_004_2_002_2  Ha Hb match to different parental hp
FM2_003_2_004_2_002_2<-data_for_selection %>%
  filter(unique_hp_FP == 2,unique_hp_3==2,unique_hp_4==2,e_match_Ha==2,e_match_Hb==2,(perc_a_m_3_Ha==max_matchperc_Ha&perc_b_m_3_Hb==max_matchperc_Hb)|(perc_a_m_3_Hb==max_matchperc_Ha&perc_b_m_3_Ha==max_matchperc_Hb))

# FM2_003_21_004_21

#FM2_003_21_004_21_002_1  Ha Hb match to same parental hp
FM2_003_21_004_21_002_1<-data_for_selection %>%
  filter(unique_hp_FP == 2,(unique_hp_3==2&unique_hp_4==1)|(unique_hp_3==1&unique_hp_4==2),
         e_match_Ha==3&e_match_Hb==3)
#FM2_003_21_004_21_002_2  Ha Hb match to different parental hp
FM2_003_21_004_21_002_2<-data_for_selection %>%
  filter(unique_hp_FP == 2,(unique_hp_3==2&unique_hp_4==1)|(unique_hp_3==1&unique_hp_4==2),
        (e_match_Ha==3&e_match_Hb==1)|(e_match_Ha==1&e_match_Hb==3))

## 2. unique_hp_FP 3

#FM3_003_21_004_21

#FM3_003_21_004_21_002_2  Ha Hb match to different parental hp
FM3_003_21_004_21_002_2<-data_for_selection %>%
  filter(unique_hp_FP == 3,(unique_hp_3==2&unique_hp_4==1)|(unique_hp_3==1&unique_hp_4==2),
        (e_match_Ha==1&e_match_Hb==2&((rep_3_Ha==2&e_match_Hb_3==2)|(rep_4_Ha==2&e_match_Hb_4==2)))|(e_match_Ha==2&e_match_Hb==1&((rep_3_Ha==2&e_match_Ha_3==2)|(rep_4_Ha==2&e_match_Ha_4==2))))


# FM3_003_2_004_2

#FM3_003_2_004_2_002_1   Ha Hb match to same parental hp
FM3_003_2_004_2_002_1<-data_for_selection %>%
  filter(unique_hp_FP == 3,unique_hp_3==2,unique_hp_4==2,
         (rep_3_Ha==2&perc_a_m_3_Ha==max_matchperc_Ha&perc_b_m_3_Ha==max_matchperc_Hb)|(rep_3_Hb==2&perc_a_m_3_Hb==max_matchperc_Ha&perc_b_m_3_Hb==max_matchperc_Hb))
#FM3_003_2_004_2_002_2_rep_1  Ha Hb match to different parental hp, they all rep 1 in parental hp  
FM3_003_2_004_2_002_2_rep_1<-data_for_selection %>%
  filter(unique_hp_FP == 3,unique_hp_3==2,unique_hp_4==2,
         e_match_Ha==1&e_match_Hb==1,
         e_match_Ha_3+e_match_Hb_3==1)                                                       
#FM3_003_2_004_2_002_2_rep_1_2  Ha Hb match to different parental hp, one rep 1 & another rep 2 in parental hp  
FM3_003_2_004_2_002_2_rep_1_2<-data_for_selection %>%
  filter(unique_hp_FP == 3,unique_hp_3==2,unique_hp_4==2,
        ((e_match_Ha==2&e_match_Hb==1)&((rep_3_Ha==2&perc_a_m_3_Ha==max_matchperc_Ha)|(rep_3_Hb==2&perc_a_m_3_Hb==max_matchperc_Ha)))|    
        ((e_match_Ha==1&e_match_Hb==2)&((rep_3_Ha==2&perc_b_m_3_Ha==max_matchperc_Hb)|(rep_3_Hb==2&perc_b_m_3_Hb==max_matchperc_Hb)))) 

# 3. unique_hp_FP 4  

#FM4_003_2_004_2_002_2  Ha Hb match to different parental hp
FM4_003_2_004_2_002_2<-data_for_selection %>%
  filter(unique_hp_FP == 4,unique_hp_3==2,unique_hp_4==2,e_match_Ha==1,e_match_Hb==1,e_match_Ha_3+e_match_Hb_3==1)

# combined successful blocks 

ncol<-ncol(FM3_003_2_004_2_002_1)+1
colnames<-colnames(FM2_003_2_004_2_002_1)

trio_success<-data.frame(matrix(0,0,ncol))
colnames(trio_success)<-c(colnames,"class")

for (category in c("FM2_003_2_004_2_002_1","FM2_003_2_004_2_002_2","FM2_003_21_004_21_002_1","FM2_003_21_004_21_002_2","FM3_003_21_004_21_002_2","FM3_003_2_004_2_002_1","FM3_003_2_004_2_002_2_rep_1","FM3_003_2_004_2_002_2_rep_1_2","FM4_003_2_004_2_002_2")){
  sub<-get(category)
  if (nrow(sub)!=0){
    sub$class<-category
    trio_success<-rbind(trio_success,sub)
  }
}

trio_success$CHROM<-as.character(trio_success$CHROM)

sorted_trio_success<-trio_success%>%
  arrange(CHROM,start)

# get HG002 Ha and Hb matched parental haplotypes respectively--> match_a match_b
trio_success_with_match_ab<-sorted_trio_success%>%
  rowwise()%>%
  mutate(match_a = paste(c("3_Ha","3_Hb","4_Ha","4_Hb")[c_across(perc_a_m_3_Ha:perc_a_m_4_Hb) == max_matchperc_Ha],collapse=",")) %>%
  mutate(match_b = paste(c("3_Ha","3_Hb","4_Ha","4_Hb")[c_across(perc_b_m_3_Ha:perc_b_m_4_Hb) == max_matchperc_Hb],collapse=","))

write.table(trio_success_with_match_ab,file=paste("b_",sample,"_trio_success.txt",sep=""),row.names = F,quote = F)

# # --------------------plotting-------------------

# select only 2 columns: len_bp,len_snv
plot_trio_success<- trio_success %>%
  select(len_bp,len_snv) 

# plot histgoram for bp length, bin=5000pb
bin_bp=5000
ggplot(plot_trio_success) +
  geom_histogram(aes(x=len_bp),binwidth = bin_bp,boundary=0) +
  scale_x_continuous(name="len_kb",labels = scales::label_number(suffix = "", scale = 1e-3))+
  theme_minimal()+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=14,face="bold"))

ggsave(
  paste(sample,"_bplen_bin_",bin_bp,".pdf",sep=""),
  plot = last_plot(),
  width = 20,
  height = 10,
  units = c("cm"),
  device = "pdf",
  bg='#ffffff')

# plot histgoram for snv length, bin=5snv
bin_snv=5
ggplot(plot_trio_success) +
  geom_histogram(aes(x=len_snv),binwidth = bin_snv,boundary=0) +
  theme_minimal()+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=14,face="bold"))


ggsave(
  paste(sample,"_SNVlen_bin_",bin_snv,".pdf",sep=""),
  plot = last_plot(),
  width = 20,
  height = 10,
  units = c("cm"),
  device = "pdf",
  bg='#ffffff')

# plot histgoram for snv length, bin=5snv, xlim=c(0,250)
ggplot(plot_trio_success) +
  geom_histogram(aes(x=len_snv),binwidth = bin_snv,boundary=0) +
  coord_cartesian(xlim = c(0,250))+
  theme_minimal()+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=14,face="bold"))


ggsave(
  paste(sample,"_SNVlen_bin_",bin_snv,"_xlim_0_250.pdf",sep=""),
  plot = last_plot(),
  width = 20,
  height = 10,
  units = c("cm"),
  device = "pdf",
  bg='#ffffff')

