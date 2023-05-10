# sum protein coding exon: 502253
# sum protein coding gene: 20028
args = commandArgs(trailingOnly=TRUE)
sample <-args[1]
setwd(args[2])
exon<-read.table(paste(sample,".exon.cover.bed",sep=""))
#minimum, lower-hinge, median, upper-hinge, maximum
stat_exon<-fivenum(exon[,3]-exon[,2])
gene<-read.table(paste(sample,".gene.cover.bed",sep=""))
stat_gene<-fivenum(gene[,3]-gene[,2])

stat<-data.frame(matrix(nrow=14,ncol=1))
rownames(stat)<-c("num_covered_pro_exon","%covered_pro_exon", "min_len_covered_pro_exon", "Q1_len_covered_pro_exon", "Q2_len_covered_pro_exon", "Q3_len_covered_pro_exon","max_len_covered_pro_exon","num_covered_pro_gene","%covered_pro_gene", "min_len_covered_pro_gene", "Q1_len_covered_pro_gene", "Q2_len_covered_pro_gene", "Q3_len_covered_pro_gene","max_len_covered_pro_gene")
stat["num_covered_pro_exon",]<-nrow(exon)
stat["%covered_pro_exon",]<-nrow(exon)/502253
stat["min_len_covered_pro_exon",]<-stat_exon[1]
stat["Q1_len_covered_pro_exon",]<-stat_exon[2]
stat["Q2_len_covered_pro_exon",]<-stat_exon[3]
stat["Q3_len_covered_pro_exon",]<-stat_exon[4]
stat["max_len_covered_pro_exon",]<-stat_exon[5]

stat["num_covered_pro_gene",]<-nrow(gene)
stat["%covered_pro_gene",]<-nrow(gene)/20028
stat["min_len_covered_pro_gene",]<-stat_gene[1]
stat["Q1_len_covered_pro_gene",]<-stat_gene[2]
stat["Q2_len_covered_pro_gene",]<-stat_gene[3]
stat["Q3_len_covered_pro_gene",]<-stat_gene[4]
stat["max_len_covered_pro_gene",]<-stat_gene[5]

write.table(stat,file=paste(sample,"_coverd_pro_exon_gene_stat.txt",sep=""),quote=F,sep="\t",row.names=T,col.names=F)

