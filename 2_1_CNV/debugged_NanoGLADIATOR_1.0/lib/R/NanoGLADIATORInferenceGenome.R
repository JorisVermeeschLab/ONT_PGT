vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


##  Setting input paths for normalized read count and experimental design ###
DataFolder <- split.vars[1]
TargetFolder <- split.vars[2]
ExperimentalFile <- split.vars[3]
ExperimentalDesign <- split.vars[4]
TargetName <- split.vars[5]
ProgramFolder <- split.vars[6]
Assembly <- split.vars[7]

options(scipen = 999)

### Load and set experimental design ###
ExperimentalTable <- read.table(ExperimentalFile,sep=" ",quote="",header=F)
LabelName <- as.character(ExperimentalTable[,1])
PathInVec <- as.character(ExperimentalTable[,2])
ExpName <- as.character(ExperimentalTable[,3])


### Create the vector for the experimental design ###
if (ExperimentalDesign=="pooling")
{
  ExpTest <- ExpName
  PathInVecTS <- PathInVec
  FileRef <- file.path(DataFolder,"Control","RCNorm","Control.NRC.RData")
  load(file=FileRef)
  DataSeqRef <- as.numeric(MatrixNorm[,6])
}

if (ExperimentalDesign=="paired")
{
  indC <- grep("C",LabelName)
  indT <- grep("T",LabelName)
  LabelNameC <- LabelName[indC]
  LabelNameT <- LabelName[indT]
  ExpNameC <- ExpName[indC]
  ExpNameT <- ExpName[indT]
  PathInVecC <- PathInVec[indC]
  PathInVecT <- PathInVec[indT]
  NumC <- as.numeric(substr(LabelNameC, 2, 100000))
  NumT <- as.numeric(substr(LabelNameT, 2, 100000))
  indCS <- sort(NumC,index.return=T)$ix
  indTS <- sort(NumT,index.return=T)$ix
  ExpTest <- ExpNameT[indTS]
  ExpControl <- ExpNameC[indCS]
  PathInVecTS <- PathInVecT[indTS]
  PathInVecCS <- PathInVecC[indCS]
}

if (ExperimentalDesign=="nocontrol")
{
  ExpTest <- c(ExpName)
  PathInVecTS <- PathInVec
}



### Loading target chromosomes ###
TargetChrom <- file.path(TargetFolder,paste(TargetName,"_chromosome.txt",sep=""))
CHR <- scan(TargetChrom,what="character")
unique.chrom <- paste(CHR,".",sep="")
Chr2Plot<-paste(CHR,collapse=",")


PythonPlot<-file.path(ProgramFolder,"/lib/python/NanoGLADIATORPlot.py")
source(file.path(ProgramFolder,"/lib/R/LibraryFractionCall.R"))
dyn.load(file.path(ProgramFolder,"/lib/F77/FastJointSLMLibrary.so"))
source(file.path(ProgramFolder,"/lib/R/LibraryJSLM.R"))

### Loading centromere file ### 
CentromereTable <- read.table(file.path(ProgramFolder,paste("data/centromere/CentromerePosition_",Assembly,".txt",sep="")),sep="\t",quote="",header=T)
CentroChr <- as.character(CentromereTable[,1])
if (nchar(CHR[1])<4)
{
  CentroChr<-substr(CentroChr, 4, 100000)
}

CentroStart <- as.numeric(CentromereTable[,2])
CentroEnd <- as.numeric(CentromereTable[,3])




for (zz in 1:length(ExpTest))
{
  ExpLabelOut<-ExpTest[zz]
  ### Loading normalized read count for reference sample when ExperimentalDesign = somatic or pooling ####
  if (ExperimentalDesign=="paired")
  {
    FileNameRef<-list.files(file.path(PathInVecCS[zz],"RCNorm"))
    FileRef <- file.path(PathInVecCS[zz],"RCNorm",FileNameRef)
    load(file=FileRef)
    RefMatrixNorm <- MatrixNorm
    DataSeqRef <- as.numeric(RefMatrixNorm[,6])
    Median1CopyRef<-round(median(DataSeqRef)/2)
  }
  
  ### Loading normalized read count for test sample ###
  FileNameTest<-list.files(file.path(PathInVecTS[zz],"RCNorm"))
  FileTest <- file.path(PathInVecTS[zz],"RCNorm",FileNameTest)
  
  #FileTest <- file.path(PathInVecTS[zz],"RCNorm",paste(ExpTest[zz],".NRC.RData",sep=""))
  
  load(file=FileTest)
  
  TestMatrixNorm <- MatrixNorm
  DataSeqTest <- as.numeric(TestMatrixNorm[,6])
  Position <- as.integer(TestMatrixNorm[,2])
  chrom <- as.character(TestMatrixNorm[,1])
  start <- as.numeric(TestMatrixNorm[,3])
  end <- as.numeric(TestMatrixNorm[,4])
  Gene <- as.character(TestMatrixNorm[,5])
  Median1CopyTest<-round(median(DataSeqTest)/2)
  
  ### Calculating Log2-ratio and lowess normalization ###
  if (ExperimentalDesign!="nocontrol")
  {
    A <- 0.5*log2(DataSeqTest*DataSeqRef)
    M <- log2(DataSeqTest/DataSeqRef)
    LogDataNorm <- M - median(M)
    
  }
  
  if (ExperimentalDesign=="nocontrol")
  {
    LogDataNorm<-log2(DataSeqTest/median(DataSeqTest))
  }
  
  
  
  
  ### Setting starting parameters of HSLM ### 
  ParameterVec <- as.numeric(unlist(read.table(file.path(ProgramFolder,"/ParameterFile.txt"),sep="\t",quote="",header=F,comment.char="#")))
  omega <- ParameterVec[1]
  eta <- ParameterVec[2]
  stepeta <- ParameterVec[3]
  cell <- ParameterVec[4]
  thrd <- ParameterVec[5]
  thru <- ParameterVec[6]
  FW <- ParameterVec[7]
  
  ### Calculating parameters of the HSLM ###
  mw <- 1
  ParamList <- ParamEstSeq(rbind(LogDataNorm),omega)
  mi <- ParamList$mi
  smu <- ParamList$smu
  sepsilon <- ParamList$sepsilon
  muk <- MukEst(rbind(LogDataNorm),mw)
  
  
  ###  Segmentation of the log2-ratio profiles with HSLM ###
  MatrixCall<-c()
  MatrixSeg <- matrix(NA,nrow=length(LogDataNorm),ncol=7)
  indCountS<-1
  for (i in 1:length(CHR))
  {
    MatrixCallChr<-c()
    chr <- as.character(unlist(CHR))[i]
    
    indchr <- which(chrom==chr)
    seqChrom <- LogDataNorm[indchr]
    PosChrom <- Position[indchr]
    startChrom <- start[indchr]
    endChrom <- end[indchr]
    GeneChrom <- Gene[indchr]
    DataSeqTestChrom <- DataSeqTest[indchr]
    if (ExperimentalDesign=="paired")
    {
      DataSeqRefChrom<-DataSeqRef[indchr]
    }
    
    indCountE<-indCountS+length(indchr)-1
    
    
    splitchrom1 <- CentroStart[which(CentroChr==chr)]
    splitchrom2 <- CentroEnd[which(CentroChr==chr)]
    
    splitind1 <- tail(which(PosChrom<splitchrom1),1)
    splitind2 <- which(PosChrom>splitchrom2)[1]
    ChrSeg<-c()
    if (length(splitind1)!=0)
    {
      ind1 <- c(1:splitind1)
      ind2 <- c(splitind2:length(seqChrom))
      
      DataSeqTestChrom1<-DataSeqTestChrom[ind1]
      DataSeqTestChrom2<-DataSeqTestChrom[ind2]
      
      if (ExperimentalDesign=="paired")
      {
        DataSeqRefChrom1<-DataSeqRefChrom[ind1]
        DataSeqRefChrom2<-DataSeqRefChrom[ind2]
      }
      
      
      
      DataSeq1 <- rbind(seqChrom[ind1])
      DataSeq2 <- rbind(seqChrom[ind2])
      Pos1 <- PosChrom[1:splitind1]
      Pos2 <- PosChrom[splitind2:length(seqChrom)]
      start1 <- startChrom[1:splitind1]
      end1 <- endChrom[1:splitind1]
      start2 <- startChrom[splitind2:length(seqChrom)]
      end2 <- endChrom[splitind2:length(seqChrom)]
      Gene1 <- GeneChrom[ind1]
      Gene2 <- GeneChrom[ind2]
      
      
      TotalPredBreak1 <- JointSeg(DataSeq1,eta,omega,muk,mi,smu,sepsilon)
      if (length(which(is.na(TotalPredBreak1)))==0)
      {
        TotalPredBreak1 <- FilterSeg(TotalPredBreak1,FW)
        if (ExperimentalDesign=="paired")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallTwoSample(TotalPredBreak1,DataSeqTestChrom1,DataSeqRefChrom1,DataSeq1,Median1CopyTest,Median1CopyRef,Pos1,chr))
          
        }
        if (ExperimentalDesign=="nocontrol")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallOneSample(TotalPredBreak1,DataSeqTestChrom1,DataSeq1,Median1CopyTest,Pos1,chr))
          
        }
        DataSeg1 <- SegResults(DataSeq1,TotalPredBreak1)
        ChrSeg<-rbind(ChrSeg,cbind(rep(chr,length(ind1)),Pos1,Gene1,start1,end1,t(DataSeq1),t(DataSeg1)))
      }
      if (length(which(is.na(TotalPredBreak1)))!=0)
      {
        print(paste("The HSLM analysis of short arm of chromosome",chr,"was aborted because the total number of EMRC is too small:",length(ind1)))
        ChrSeg<-rbind(ChrSeg,cbind(rep(chr,length(ind1)),Pos1,Gene1,start1,end1,matrix(NA,nrow=length(ind1),ncol=2)))
      }
      #####-------------add also centromere region in ChrSeg-----------
      indCountcentro<-which(PosChrom>=splitchrom1 & PosChrom<=splitchrom2)
      Poscen <- PosChrom[indCountcentro]
      startcen <- startChrom[indCountcentro]
      endcen <- endChrom[indCountcentro]
      Genecen <- GeneChrom[indCountcentro]

      ChrSeg<-rbind(ChrSeg,cbind(rep(chr,length(indCountcentro)),Poscen,Genecen,startcen,endcen,matrix(NA,nrow=length(indCountcentro),ncol=2)))
     ####---------------------------------------------------------------

      TotalPredBreak2 <- JointSeg(DataSeq2,eta,omega,muk,mi,smu,sepsilon)
      if (length(which(is.na(TotalPredBreak1)))==0)
      {
        TotalPredBreak2 <- FilterSeg(TotalPredBreak2,FW)
        if (ExperimentalDesign=="paired")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallTwoSample(TotalPredBreak2,DataSeqTestChrom2,DataSeqRefChrom2,DataSeq2,Median1CopyTest,Median1CopyRef,Pos2,chr))

        }
        if (ExperimentalDesign=="nocontrol")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallOneSample(TotalPredBreak2,DataSeqTestChrom2,DataSeq2,Median1CopyTest,Pos2,chr))
          
        }
        DataSeg2 <- SegResults(DataSeq2,TotalPredBreak2)
        ChrSeg<-rbind(ChrSeg,cbind(rep(chr,length(ind2)),Pos2,Gene2,start2,end2,t(DataSeq2),t(DataSeg2)))
      }
      if (length(which(is.na(TotalPredBreak2)))!=0)
      {
        print(paste("The HSLM analysis of long arm of chromosome",chr,"was aborted because the total number of EMRC is too small:",length(ind1)))
        ChrSeg<-rbind(ChrSeg,cbind(rep(chr,length(ind2)),Pos2,Gene2,start2,end2,matrix(NA,nrow=length(ind2),ncol=2)))
      }
    }
    
    if (length(splitind1)==0)
    {
      seqChrom <- rbind(seqChrom)
      TotalPredBreak <- JointSeg(seqChrom,eta,omega,muk,mi,smu,sepsilon)
      if (length(which(is.na(TotalPredBreak)))==0)
      {
        TotalPredBreak <- FilterSeg(TotalPredBreak,FW)
        if (ExperimentalDesign=="paired")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallTwoSample(TotalPredBreak,DataSeqTestChrom,DataSeqRefChrom,seqChrom,Median1CopyTest,Median1CopyRef,PosChrom,chr))
        }
        
        if (ExperimentalDesign=="nocontrol")
        {
          MatrixCallChr<-rbind(MatrixCallChr,FractionCallOneSample(TotalPredBreak,DataSeqTestChrom,seqChrom,Median1CopyTest,PosChrom,chr))
          
        }
        DataSeg <- SegResults(seqChrom,TotalPredBreak)
        ChrSeg <- cbind(rep(chr,length(indchr)),PosChrom,cbind(GeneChrom),startChrom,endChrom,t(seqChrom),t(DataSeg))
      }
      if (length(which(is.na(TotalPredBreak)))!=0)
      {
        print(paste("The HSLM analysis of chromosome",chr,"was aborted because the total number of EMRC is too small:",length(indchr)))
        ChrSeg <- rbind(ChrSeg,cbind(rep(chr,length(indchr)),PosChrom,cbind(GeneChrom),startChrom,endChrom,matrix(NA,nrow=length(indchr),ncol=2)))
      }
    }
    
    ##### Saving Segmentation Chrom by Chrom ###
    ChrSegOut <- ChrSeg[,c(1,2,4,5,6,7)]
    ChrSegOut <- na.omit(ChrSegOut)
    colnames(ChrSegOut) <- c("Chromosome","Position","Start","End","Log2R","SegMean")
    FileOutSegChr <- file.path(DataFolder,"Results",ExpTest[zz],paste(chr,"_HSLMResults_",ExpTest[zz],".txt",sep=""))
    write.table(ChrSegOut,FileOutSegChr,col.names=T,row.names=F,sep="\t",quote=F)
    
    ##### Saving Fraction Calling Chrom by Chrom ###
    FileOutCallChr <- file.path(DataFolder,"Results",ExpTest[zz],paste(chr,"_FractionCallResults_",ExpTest[zz],".txt",sep=""))
    write.table(MatrixCallChr,FileOutCallChr,col.names=T,row.names=F,sep="\t",quote=F)
    
    
    
    MatrixCall<-rbind(MatrixCall,MatrixCallChr)
    MatrixSeg[c(indCountS:indCountE),] <- ChrSeg
    indCountS<-indCountE+1
  }
  
  
  #### Saving Segmented Matrix ####
  MatrixSeg1 <- MatrixSeg[,c(1,2,4,5,6,7)]
  colnames(MatrixSeg1) <- c("Chromosome","Position","Start","End","Log2R","SegMean")
  FileOutSeg <- file.path(DataFolder,"Results",ExpTest[zz],paste("HSLMResults_",ExpTest[zz],".txt",sep=""))
  write.table(MatrixSeg1,FileOutSeg,col.names=T,row.names=F,sep="\t",quote=F)
    
    
 
  
  FileOutCall <- file.path(DataFolder,"Results",ExpTest[zz],paste("FractionCallResults_",ExpTest[zz],".txt",sep=""))
  write.table(MatrixCall,FileOutCall,col.names=F,row.names=F,sep="\t",quote=F)
  
  ##### Plotting Results in html ####
  FilePathPlotIn<-file.path(DataFolder,"Results",ExpTest[zz])
  FilePathPlotOut<-file.path(DataFolder,"Plots",ExpTest[zz])
  # Python2Run<-paste("python ",PythonPlot," ",ExpLabelOut," ",Chr2Plot," ",FilePathPlotIn," ",FilePathPlotOut,sep="")
    Python2Run<-paste("/vsc-hard-mounts/leuven-data/339/vsc33900/miniconda3/bin/python ",PythonPlot," ",ExpLabelOut," ",Chr2Plot," ",FilePathPlotIn," ",FilePathPlotOut,sep="")

  system(Python2Run,intern=FALSE)
  
  
  #### FastCall Results in VCF format for Regions and Windows #####
  #VCFRegionCreate(Assembly,DataFolder,ExpLabelOut,TargetFolder,SummaryData,MetaData,out)
}





