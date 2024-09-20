library(dplyr)
library(ggplot2)
library(PSCBS)
library(ggpubr)
library(tidyr)

sample<-""
name<-""
setwd("")

# -----get the length of each chr-----------
ChrsLength<-read.table("xx/hg38.bed")
colnames(ChrsLength)<-c("chr","start","end")
ChrsLength <- ChrsLength%>%
  mutate(seq=case_when(chr=="X"~"23",
                       TRUE~chr))%>%
  arrange(as.numeric(seq))%>%
  select(chr,end)

CumLengths <- ChrsLength
CumLengths[,"Length"] <- cumsum(as.numeric(ChrsLength[,2]))
GenomeLength <- as.numeric(CumLengths[CumLengths[,"chr"]=="X","Length"])

# ----------read in and prepare data-------

data<-read.csv("diff_ref_vcf_info.txt",header = T,sep="\t")%>%
  filter(DP_E>=5&DP_E<=50)%>%
  select(-DP_E,-DP_F,-DP_M)

# change AF 
AF<-data%>%
  mutate(plot_AF_Father_E=case_when(GT_F==0~1-AF_E,
         TRUE~AF_E))%>%
  mutate(plot_AF_Mother_E=case_when(GT_M==0~1-AF_E,
                                      TRUE~AF_E))%>%
  select("CHROM","POS","plot_AF_Father_E","plot_AF_Mother_E")

AF_change_pos<-AF
chrs<-unique(AF$CHROM)

for(chr in chrs){
  if (chr != 1){
    ToAdd <- as.numeric(CumLengths[grep(paste("^",chr,"$",sep=""),CumLengths[,"chr"])-1,"Length"])
    AF_change_pos[AF_change_pos$CHROM==chr,"POS"] <- AF[AF$CHROM==chr,"POS"] +ToAdd
  }
}
# ------CBS segmentation ; binsize 1MB ------

bin_means<-AF_change_pos%>%  
  group_by(CHROM)%>% 
  mutate(bin_start=1000000*((POS-1) %/% 1000000)+1)%>%
  mutate(bin_end=1000000*((POS-1) %/% 1000000+1))%>%
  group_by(bin_start) %>%
  summarise(CHROM=unique(CHROM),
            bin_start=unique(bin_start),
            bin_end=unique(bin_end),
            mean_AF_Father_E=mean(plot_AF_Father_E),
            mean_AF_Mother_E=mean(plot_AF_Mother_E),
            dens_loci=1000000/n()) %>%
  distinct()%>%
  filter(dens_loci<=30000)

chrs=unique(CumLengths$chr)
bin_means_revised<-bin_means
for (chr in chrs){
  chr_end=CumLengths[CumLengths$chr==chr,"Length"]
  chr_start=chr_end-CumLengths[chr,"end"]+1
  bin_means_revised<-bin_means_revised%>%
    mutate(bin_start=case_when(CHROM==chr&bin_start<chr_start~chr_start,
                                TRUE~bin_start))%>%
             mutate(bin_end=case_when(CHROM==chr&bin_end>chr_end~chr_end,
                                         TRUE~bin_end))
                               
}

bin_means_revised<-bin_means_revised%>%
  mutate(bin=paste0(CHROM,":",bin_start,"-",bin_end))

# mean_AF_Father_E segmentation
data <- bin_means_revised[, c("CHROM", "bin_start","mean_AF_Father_E")]
data$CHROM[data$CHROM=="X"]<-23
colnames(data) <-c("chromosome","x", "y")
# drops single-locus outliers along the genome that have a signal that differ significantly from the neighboring loci.
data <- dropSegmentationOutliers(data)
#Identifying TCN segments
fit <- segmentByCBS(data)
Segments_F<- getSegments(fit, simplify = TRUE)%>%
  select(chromosome,start,end,mean)%>%
  drop_na()%>%
  mutate(start=start,end=end+1000000-1)

Segments_F$chromosome[Segments_F$chromosome=="23"]<-"X"

for (chr in chrs){
  chr_end=CumLengths[CumLengths$chr==chr,"Length"]
  chr_start=chr_end-CumLengths[chr,"end"]+1
  Segments_F<-Segments_F%>%
    mutate(start=case_when(chromosome==chr&start<chr_start~chr_start,
                           TRUE~start))%>%
    mutate(end=case_when(chromosome==chr&end>chr_end~chr_end,
                         TRUE~end))
  
}


# mean_AF_Mother_E segmentation

data <- bin_means_revised[, c("CHROM", "bin_start","mean_AF_Mother_E")]
data$CHROM[data$CHROM=="X"]<-23
colnames(data) <-c("chromosome","x", "y")
# drops single-locus outliers along the genome that have a signal that differ significantly from the neighboring loci.
data <- dropSegmentationOutliers(data)
#Identifying TCN segments
fit <- segmentByCBS(data)
Segments_M<- getSegments(fit, simplify = TRUE)%>%
  select(chromosome,start,end,mean)%>%
  drop_na()%>%
  mutate(start=start,end=end+1000000-1)

Segments_M$chromosome[Segments_M$chromosome=="23"]<-"X"

for (chr in chrs){
  chr_end=CumLengths[CumLengths$chr==chr,"Length"]
  chr_start=chr_end-CumLengths[chr,"end"]+1
  Segments_M<-Segments_M%>%
    mutate(start=case_when(chromosome==chr&start<chr_start~chr_start,
                           TRUE~start))%>%
    mutate(end=case_when(chromosome==chr&end>chr_end~chr_end,
                         TRUE~end))
  
}

# write segmentation to file 
colnames(Segments_F)[4]<-"seg"
colnames(Segments_M)[4]<-"seg"
write.table(Segments_F,"Segments_F.txt",row.names = F,quote = F)
write.table(Segments_M,"Segments_M.txt",row.names = F,quote = F)


# plot 
line_width_v=0.1
line_width_h=0.2
line_width_seg=0.5
p_M<-ggplot()+
  geom_point(data=bin_means_revised,aes(x=(bin_start+bin_end)/2,y=mean_AF_Mother_E),colour="grey",size=1,alpha=0.1)+
  geom_segment(data=Segments_M,aes(x=start,xend=end,y=seg,yend=seg),color="red",linewidth=line_width_seg)+
  geom_vline(data=CumLengths,aes(xintercept = Length),linetype="dotted",color="black",linewidth=line_width_v)+
  geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=line_width_v)+
  geom_hline(yintercept = 0.33,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.5,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.67,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 1,linetype="dotted",color="black",linewidth=line_width_h)+
  labs(title=name,
       x="",
       y="Mat AF")+
  scale_y_continuous(limits =c(-0.1,1), breaks = c(0,0.5,1))+
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0,GenomeLength), 
                     breaks = CumLengths$Length-CumLengths$end/2, 
                     labels = CumLengths$chr)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title=element_text(size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(r = 0)),
        axis.text=element_text(size=4,face="bold"),
        axis.title=element_text(size=8,face="bold"))

p_F<-ggplot()+
  geom_point(data=bin_means_revised,aes(x=(bin_start+bin_end)/2,y=mean_AF_Father_E),colour="grey",size=1,alpha=0.1)+
  geom_segment(data=Segments_F,aes(x=start,xend=end,y=seg,yend=seg),colour="blue",linewidth=line_width_seg)+
  geom_vline(data=CumLengths,aes(xintercept = Length),linetype="dotted",color="black",linewidth=line_width_v)+
  geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=line_width_v)+
  geom_hline(yintercept = 0.33,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.5,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.67,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 1,linetype="dotted",color="black",linewidth=line_width_h)+
  labs(title="",
       x="Chromosome",
       y="Pat AF")+
  scale_y_continuous(limits =c(-0.1,1), breaks = c(0,0.5,1))+
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0,GenomeLength), 
                     breaks = CumLengths$Length-CumLengths$end/2, 
                     labels = CumLengths$chr)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title=element_text(size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(r = 0)),
        axis.text=element_text(size=4,face="bold"),
        axis.title=element_text(size=8,face="bold"))

figure <- ggarrange(p_M, p_F,
                    ncol =1, nrow = 2)

ggsave(
  paste0(sample,"_parental_AF.jpg"),
  plot = last_plot(), 
  device = "jpg",
  path = NULL,
  scale = 1,
  width = 150,
  height = 80,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
)

  