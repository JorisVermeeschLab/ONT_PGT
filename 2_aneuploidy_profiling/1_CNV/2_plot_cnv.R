library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
sample=args[2]

# -----get the length of each chr-----------
ChrsLength<-read.table("/Users/miffy/OneDrive - KU Leuven/1_working directory/02_nanopore/00_data_for_manuscript/0_plotting_for_manuscript/3_2_CNV_plots_hg38w1000000/run_plot_R/hg38.bed")
colnames(ChrsLength)<-c("chr","start","end")
ChrsLength <- ChrsLength%>%
  mutate(seq=case_when(chr=="X"~"23",
                       TRUE~chr))%>%
  arrange(as.numeric(seq))%>%
  select(chr,end)

CumLengths <- ChrsLength
CumLengths[,"Length"] <- cumsum(as.numeric(ChrsLength[,2]))
GenomeLength <- as.numeric(CumLengths[CumLengths[,"chr"]=="X","Length"])

# -----read infercnv results -----------

cnv<-read.table(args[3],header = T)%>%
  filter(Chromosome!="Y")

  cnv_change_pos<-cnv
  chrs<-unique(cnv$Chromosome)[1:23]
  
  for(chr in chrs){
    if (chr != 1){
      ToAdd <- as.numeric(CumLengths[grep(paste("^",chr,"$",sep=""),CumLengths[,"chr"])-1,"Length"])
      cnv_change_pos[cnv_change_pos$Chromosome==chr,"Start"] <- cnv[cnv$Chromosome==chr,"Start"] +ToAdd
      cnv_change_pos[cnv_change_pos$Chromosome==chr,"End"] <- cnv[cnv$Chromosome==chr,"End"] +ToAdd
      cnv_change_pos[cnv_change_pos$Chromosome==chr,"Position"] <- cnv[cnv$Chromosome==chr,"Position"] +ToAdd
    }
  }

seg<-cnv_change_pos%>%
  select(Chromosome,Start,End,SegMean)%>%
  distinct()
line_width_v=0.1
line_width_h=0.2
line_width_seg=0.5
ggplot()+
  geom_point(data=cnv_change_pos,aes(x=Position,y=Log2R),colour="grey",size=0.3,alpha=0.2)+
  geom_segment(data=seg,aes(x=Start,xend=End,y=SegMean,yend=SegMean),colour="orange",linewidth=line_width_seg)+
  scale_fill_identity()+
  geom_vline(data=CumLengths,aes(xintercept = Length),linetype="dotted",color="black",linewidth=line_width_v)+
  geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=line_width_v)+
  geom_hline(yintercept = log2(3/2),linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = log2(2/2),linetype="dotted",color="black",linewidth=line_width_h)+
  geom_hline(yintercept = log2(1/2),linetype="dotted",color="black",linewidth=line_width_h)+
  labs(title=args[4],
     x="Chromosome",
     y="Log2 ratio")+
  ylim(c(-1.8,1.8))+
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

ggsave(
  paste0(sample,"_cnv.jpg"),
  plot = last_plot(), 
  device = "jpg",
  path = NULL,
  scale = 1,
  width = 150,
  height = 40,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
)

  
  
  
  
  