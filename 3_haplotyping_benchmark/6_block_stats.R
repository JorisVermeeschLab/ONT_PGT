library(dplyr)
library(ggplot2)
library(tidyr)


args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

file<-paste0("b_",sample,"_label.txt")

  name_0<-sub("_label\\.txt","",file)
  name<-sub("b_","",name_0)
  data_raw<-read.table(file,header = T)
  data<-data_raw%>%
    select("ID","len_bp","len_snv")
  data$name<-name

# Length (kb)
ggplot(data,aes(x = len_bp/1000)) + 
  geom_histogram(binwidth = 20)+
  facet_grid(~name)+
  labs(y="Count",
       x="Length (kb)")+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.text=element_text(size=16,face="bold"),
    axis.title=element_text(size=22,face="bold"),
    strip.text = element_text(size = 16, face="bold")
  )

ggsave(
  "length hist resolved unchanged blocks bin20kb.jpg",
  plot = last_plot(), 
  device = "jpg",
  path = NULL,
  scale = 1,
  width = 270,
  height = 100,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
) 

# Length iSNP
ggplot(data,aes(x = len_snv)) + 
  geom_histogram(binwidth = 20)+
  facet_grid(~name)+
  labs(y="Count",
       x="Number of iSNVs")+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.text=element_text(size=16,face="bold"),
    axis.title=element_text(size=22,face="bold"),
    strip.text = element_text(size = 16, face="bold")
  )

ggsave(
  "hist iSNV count bin20.jpg",
  plot = last_plot(), 
  device = "jpg",
  path = NULL,
  scale = 1,
  width = 270,
  height = 100,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
) 

# output summary stats

summary_stats<-data%>%
  group_by(name)%>%
  summarize(n_units=n(),
            mean_len_bp=mean(len_bp),
            median_len_bp=median(len_bp),
            sd_len_bp=sd(len_bp),
            max_len_bp=max(len_bp),
            min_len_bp=min(len_bp),
            mean_len_snv=mean(len_snv),
            median_len_snv=median(len_snv),
            sd_len_snv=sd(len_snv),
            max_len_snv=max(len_snv),
            min_len_snv=min(len_snv),
            )%>%
            t()
  
write.table(summary_stats,file="f_sumary_stats.txt",col.names = F,sep="\t",quote = F)
  
