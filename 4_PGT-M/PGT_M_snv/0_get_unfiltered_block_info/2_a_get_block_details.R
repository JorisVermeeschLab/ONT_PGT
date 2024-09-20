library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# read in 
data_raw<-read.table("parents_vcf_info_filtered_revised.txt",header = T)

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
       filter(!((CHROM_new=="PAR1"|CHROM_new=="PAR2")&(GT_HG003=="0|0"&GT_HG004=="0|0")))

# change PS on chrX specific regions for HG003 --> make the 3 regions distinct. 
# same paternal PS -- > blocks are separated by only maternal PS for chrx specific regions.

data<-data_MIC_labeled %>% 
      mutate(PS_HG003=case_when(CHROM_new=="X"&POS<10001 ~ "1",
                          CHROM_new=="X"&POS>2781479&POS<155701383 ~ "2781480",
                          CHROM_new=="X"&POS>156030895 ~ "156030896",
                          TRUE ~ PS_HG003)) 

 # replace . in PS_** columns with value from the surrounding line
data_ps_imputed<-data %>%
  group_by(CHROM_new) %>%
  mutate(PS_HG003 = as.numeric(na_if(PS_HG003, '.'))) %>%
  mutate(PS_HG004 = as.numeric(na_if(PS_HG004, '.'))) %>%
  fill(PS_HG003,PS_HG004,.direction = "downup") %>% 
  ungroup() 

# group by CHROM_new & PS; gave each group a unique ID; PS_HG003+PS_HG004 enough to label a block
# numbering to ensure length<1MB 
# seperate GT 
data_gt_sep_ps_imputed<-data_ps_imputed %>%                   
  mutate(numbering = POS %/% 1000000+1) %>%
  mutate(ID=paste(CHROM_new,PS_HG003,PS_HG004,numbering,sep="_")) %>%                      
  separate(GT_HG003,c("HG003_Ha","HG003_Hb"),"\\|")  %>%
  separate(GT_HG004,c("HG004_Ha","HG004_Hb"),"\\|") 

write.table(data_gt_sep_ps_imputed,file="a_FM_gt_sep_ps_imputed.txt",row.names = F,quote = F)
