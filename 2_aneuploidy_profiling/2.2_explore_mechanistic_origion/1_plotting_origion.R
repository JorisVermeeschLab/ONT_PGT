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
data<-read.table("check_origion.txt",header = T)%>%
  filter(CHROM!="Y")

# change pos
data_change_pos<-data
chrs<-unique(data$CHROM)

for(chr in chrs){
  if (chr != 1){
    ToAdd <- as.numeric(CumLengths[grep(paste("^",chr,"$",sep=""),CumLengths[,"chr"])-1,"Length"])
    data_change_pos[data_change_pos$CHROM==chr,"POS"] <- data[data$CHROM==chr,"POS"] +ToAdd
  }
}

# revise 
revised<-data_change_pos%>%
  filter(GT_E=="0/1")%>%
  separate(GT_E, into = c("Ex","Ey"), sep = "\\/",convert = TRUE) %>%
  separate(GT_F, into = c("Fx","Fy"), sep = "\\/",convert = TRUE) %>%
  separate(GT_M, into = c("Mx","My"), sep = "\\/",convert = TRUE) %>%
  select(-REF,-ALT)

revised$AF_E<-as.numeric(revised$AF_E)

# get fraction for pat mat HI; assume diff allele as H1
data_chr_F<-revised%>%
  filter((Fx==0&Fy==1)|(Fx==1&Fy==0)) %>%
  mutate(Pat_H1=case_when((Ex==0&Ey==0)&(Mx==0&My==0)~0,
                           ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Mx==0&My==0)~1,
                           (Ex==1&Ey==1)&(Mx==1&My==1)~1,
                           ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Mx==1&My==1)~0))%>%
  mutate(Pat_H2=case_when((Ex==0&Ey==0)&(Mx==0&My==0)~1,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Mx==0&My==0)~0,
                          (Ex==1&Ey==1)&(Mx==1&My==1)~0,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Mx==1&My==1)~1))%>%
  mutate(mat_hap=case_when(Mx==0&My==0~0,
                           Mx==1&My==1~1))%>%
  mutate(fraction_Pat_H1=case_when(Pat_H1==0~(1-AF_E),
                                    Pat_H1==1~AF_E))


data_chr_M<-revised%>%
  filter((Mx==0&My==1)|(Mx==1&My==0)) %>%
  mutate(Mat_H1=case_when((Ex==0&Ey==0)&(Fx==0&Fy==0)~0,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Fx==0&Fy==0)~1,
                          (Ex==1&Ey==1)&(Fx==1&Fy==1)~1,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Fx==1&Fy==1)~0))%>%
  mutate(Mat_H2=case_when((Ex==0&Ey==0)&(Fx==0&Fy==0)~1,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Fx==0&Fy==0)~0,
                          (Ex==1&Ey==1)&(Fx==1&Fy==1)~0,
                          ((Ex==0&Ey==1)|(Ex==1&Ey==0))&(Fx==1&Fy==1)~1))%>%
  mutate(pat_hap=case_when(Fx==0&Fy==0~0,
                           Fx==1&Fy==1~1))%>%
  mutate(fraction_Mat_H1=case_when(Mat_H1==0~(1-AF_E),
                                   Mat_H1==1~AF_E))

# ------CBS segmentation Father ; binsize 1MB ------

bin_means<-data_chr_F%>%  
  group_by(CHROM)%>% 
  mutate(bin_start=1000000*((POS-1) %/% 1000000)+1)%>%
  mutate(bin_end=1000000*((POS-1) %/% 1000000+1))%>%
  group_by(bin_start) %>%
  summarise(CHROM=unique(CHROM),
            bin_start=unique(bin_start),
            bin_end=unique(bin_end),
            mean_fraction_Pat_H1=mean(fraction_Pat_H1),
            dens_loci=1000000/n()) %>%
  distinct()%>%
  filter(dens_loci<=30000)

# cope with terminal bins
chrs=unique(CumLengths$chr)
bin_means_revised_F<-bin_means
for (chr in chrs){
  chr_end=CumLengths[CumLengths$chr==chr,"Length"]
  chr_start=chr_end-CumLengths[chr,"end"]+1
  bin_means_revised_F<-bin_means_revised_F%>%
    mutate(bin_start=case_when(CHROM==chr&bin_start<chr_start~chr_start,
                               TRUE~bin_start))%>%
    mutate(bin_end=case_when(CHROM==chr&bin_end>chr_end~chr_end,
                             TRUE~bin_end))
  
}

bin_means_revised_F<-bin_means_revised_F%>%
  mutate(bin=paste0(CHROM,":",bin_start,"-",bin_end))

# mean segmentation
data <- bin_means_revised_F[, c("CHROM", "bin_start","mean_fraction_Pat_H1")]
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

# ------CBS segmentation Mother ; binsize 1MB ------
bin_means<-data_chr_M%>%  
  group_by(CHROM)%>% 
  mutate(bin_start=1000000*((POS-1) %/% 1000000)+1)%>%
  mutate(bin_end=1000000*((POS-1) %/% 1000000+1))%>%
  group_by(bin_start) %>%
  summarise(CHROM=unique(CHROM),
            bin_start=unique(bin_start),
            bin_end=unique(bin_end),
            mean_fraction_Mat_H1=mean(fraction_Mat_H1),
            dens_loci=1000000/n()) %>%
  distinct()%>%
  filter(dens_loci<=30000)

# cope with terminal bins
chrs=unique(CumLengths$chr)
bin_means_revised_M<-bin_means
for (chr in chrs){
  chr_end=CumLengths[CumLengths$chr==chr,"Length"]
  chr_start=chr_end-CumLengths[chr,"end"]+1
  bin_means_revised_M<-bin_means_revised_M%>%
    mutate(bin_start=case_when(CHROM==chr&bin_start<chr_start~chr_start,
                               TRUE~bin_start))%>%
    mutate(bin_end=case_when(CHROM==chr&bin_end>chr_end~chr_end,
                             TRUE~bin_end))
  
}

bin_means_revised_M<-bin_means_revised_M%>%
  mutate(bin=paste0(CHROM,":",bin_start,"-",bin_end))

# mean segmentation
data <- bin_means_revised_M[, c("CHROM", "bin_start","mean_fraction_Mat_H1")]
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

# rm chr x 

plotM_1<-bin_means_revised_M%>%
  filter(CHROM!="X")
plotM_2<-Segments_M%>%
  filter(chromosome!="X")


plotF_1<-bin_means_revised_F%>%
  filter(CHROM!="X")
plotF_2<-Segments_F%>%
  filter(chromosome!="X")
# theis method can't check composition of aneuploidy  from single parents
# --> check ADO ratio in such case --> should be high if all are the same haplotype

p_M<-ggplot()+
  geom_point(data=plotM_1,aes(x=(bin_start+bin_end)/2,y=mean_fraction_Mat_H1),colour="grey",size=1,alpha=0.1)+
  geom_segment(data=plotM_2,aes(x=start,xend=end,y=seg,yend=seg),color="red",linewidth=line_width_seg)+
  geom_vline(data=CumLengths,aes(xintercept = Length),linetype="dotted",color="black",linewidth=line_width_v)+
  geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=line_width_v)+
  geom_hline(yintercept = 0.33,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.5,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.67,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 1,linetype="dotted",color="black",linewidth=line_width_h)+
  labs(title=name,
       x="",
       y="UAF")+
  scale_y_continuous(limits =c(-0.1,1), breaks = c(0,1/3,0.5,2/3,1),labels =c(0,0.33,0.5,0.67,1))+
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
  geom_point(data=plotF_1,aes(x=(bin_start+bin_end)/2,y=mean_fraction_Pat_H1),colour="grey",size=1,alpha=0.1)+
  geom_segment(data=plotF_2,aes(x=start,xend=end,y=seg,yend=seg),colour="blue",linewidth=line_width_seg)+
  geom_vline(data=CumLengths,aes(xintercept = Length),linetype="dotted",color="black",linewidth=line_width_v)+
  geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=line_width_v)+
  geom_hline(yintercept = 0.33,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.5,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0.67,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 0,linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = 1,linetype="dotted",color="black",linewidth=line_width_h)+
  labs(title="",
       x="Chromosome",
       y="Pat H1 fraction")+
  scale_y_continuous(limits =c(-0.1,1), breaks = c(0,1/3,0.5,2/3,1),labels =c(0,0.33,0.5,0.67,1))+
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
  paste0(sample,"_aneuploidy_origion.jpg"),
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

