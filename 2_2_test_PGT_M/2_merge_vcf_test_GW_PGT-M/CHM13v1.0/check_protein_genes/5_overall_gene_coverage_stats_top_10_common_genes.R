# get gene coverage and related stats for each SNP lim and  max_match_perc lim pair

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]
num_gene_ID_protein_coding=as.numeric(args[3])

file<-paste("b_",sample,"_trio_success.txt",sep="")
data<-read.table(file,header = T)

# plot all covered genes
covered_stats<-data %>%
  group_by(gene_id) %>%
  summarise(CHROM=unique(CHROM),
            gene_name=unique(gene_name),
            gene_id=unique(gene_id),
            gene_len_bp=unique(gene_len_bp),
            sum_len_bp=sum(len_bp),
            sum_len_snv_=sum(len_snv),
            perc_covered=sum(len_bp)/unique(gene_len_bp))

num_covered<-nrow(covered_stats)

perc_num_genes_covered<-round(num_covered/num_gene_ID_protein_coding,3)*100
title<-paste("total genes in ref: ", num_gene_ID_protein_coding,"; shown are covered genes (",num_covered,",",perc_num_genes_covered,"%)",sep="")

ggplot(covered_stats,aes(x=gene_len_bp,y = perc_covered)) + 
  geom_point(alpha=0.3)+
  scale_x_continuous(labels = scales::comma,limits=c(0,NA))+
  #facet_wrap(~gene_len_bp<=100000,scales = "free_x")+
  labs(title=title)+
  theme_bw()+
  theme(
    text = element_text("Helvetica"),    
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  )

ggsave(
  "perc_covered_vs_gene_len_bp.pdf",
  plot = last_plot(), 
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 300,
  height = NA,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
)

ggplot(covered_stats,aes(x=perc_covered)) + 
  geom_histogram(bins=100)+
  theme_bw()+
  theme(
    text = element_text("Helvetica"),    
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  )

ggsave(
  "perc_covered_hist.pdf",
  plot = last_plot(), 
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 300,
  height = NA,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
)

# plot top 10 PGT-M genes 
top_gene_list=c("BRCA2","BRCA1","PKD1","NF1","FMR1","CFTR","DMPK","DUX4","HBB","DMD")

top_genes<-covered_stats%>%
  filter(gene_name %in% top_gene_list)

write.table(top_genes,file="stats_top_10_pgt_m_genes.txt",quote = F,row.names = F,col.names = T )




