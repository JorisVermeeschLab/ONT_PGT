library(ggplot2)
library(dplyr)

line_size=0.3
textsize=6
size=0.3
a=1

# chr2 

sample="ONT1-E03"
disease_locus=47445655


  data<-read.table(paste0("d_variant_unit_detail_1.txt"),header=T)

  # color locus of interest by called variants in parents; locus 1 
  colors<-c("grey","red")
  # F
  color_Fh1=colors[data[data$POS==disease_locus,]$F_Ha+1]
  color_Fh2=colors[data[data$POS==disease_locus,]$F_Hb+1]
  # M
  color_Mh1=colors[data[data$POS==disease_locus,]$M_Ha+1]
  color_Mh2=colors[data[data$POS==disease_locus,]$M_Hb+1]
  
data$E_Ha<-as.character(data$E_Ha)
data$E_Hb<-as.character(data$E_Hb)
data$F_Ha<-as.character(data$F_Ha)
data$F_Hb<-as.character(data$F_Hb)
data$M_Ha<-as.character(data$M_Ha)
data$M_Hb<-as.character(data$M_Hb)

ggplot(data)+
  geom_segment(aes(x=POS,xend=POS,y=0,yend=1,color=E_Ha),linewidth=size,alpha=a)+
  geom_segment(aes(x=POS,xend=POS,y=1.2,yend=2.2,color=F_Ha),linewidth=size,alpha=a)+
  geom_segment(aes(x=POS,xend=POS,y=2.4,yend=3.4,color=F_Hb),linewidth=size,alpha=a)+
  geom_segment(aes(x=POS,xend=POS,y=4.4,yend=5.4,color=E_Hb),linewidth=size,alpha=a)+
  geom_segment(aes(x=POS,xend=POS,y=5.6,yend=6.6,color=M_Ha),linewidth=size,alpha=a)+
  geom_segment(aes(x=POS,xend=POS,y=6.8,yend=7.8,color=M_Hb),linewidth=size,alpha=a)+
  scale_color_manual(breaks=c("1","0"),values=c("black","grey"), labels = c("Alt allele","Ref allele"))+
  scale_y_continuous(breaks=c(0.5,1.7,2.9,4.9,6.1,7.3),labels=c("Embryo pat H2","Pat H1","Pat H2","Embryo mat H1","Mat H1","Mat H2"))+
  scale_x_continuous(labels = scales::comma)+
  geom_segment(aes(x=disease_locus,xend=disease_locus,y=1.2,yend=2.2),color=color_Fh1,linewidth=size,alpha=a)+
  geom_segment(aes(x=disease_locus,xend=disease_locus,y=2.4,yend=3.4),color=color_Fh2,linewidth=size,alpha=a)+
  geom_segment(aes(x=disease_locus,xend=disease_locus,y=5.6,yend=6.6),color=color_Mh1,linewidth=size,alpha=a)+
  geom_segment(aes(x=disease_locus,xend=disease_locus,y=6.8,yend=7.8),color=color_Mh2,linewidth=size,alpha=a)+
  labs(x="Position on chr2",
       y="",
       title = "ONT1-E03")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    text = element_text("Helvetica"),
    strip.background = element_blank(),
    axis.ticks.x=element_line(linewidth= line_size),
    axis.ticks.y = element_blank(),
    legend.title=element_blank(),
    legend.position="right",
    legend.key.size = unit(2, 'mm'),
    legend.box.spacing = unit(0, "mm"),# The spacing between the plotting area and the legend box (unit)
    axis.text=element_text(size=textsize),
    axis.text.y = element_text(colour = c("black","red","black","black","black","black")),
    axis.title = element_text(size = textsize),
    strip.text = element_text(size = textsize),
    plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
    legend.text=element_text(size=textsize,margin = margin(l = -2, unit = "mm"))
  )

ggsave(
  paste0(sample," PGT-M.jpg"),
  plot = last_plot(), 
  device = "jpg",
  path = NULL,
  scale = 1,
  width = 100,
  height =40,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg='#ffffff' #set the panel background to blank
)



