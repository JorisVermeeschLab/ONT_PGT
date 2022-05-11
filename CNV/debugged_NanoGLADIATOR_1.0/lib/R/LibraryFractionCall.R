FractionCallOneSample<-function(TotalPredBreak,DataSeqTestChrom,DataSeq,Median1CopyTest,Pos,chr)
{
  MatrixOut<-c()
  for (ll in 1:(length(TotalPredBreak)-1))
  {
    indSeg<-c((TotalPredBreak[ll]+1):TotalPredBreak[ll+1])
    NWindows<-length(indSeg)
    SegValue<-median(DataSeq[indSeg])
    TotalReadsSeg<-sum(DataSeqTestChrom[indSeg])
        
    if (SegValue>= 0.9)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Multi-Copy Duplication",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    if (SegValue<= -1.2)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Double Deletion",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    
    if (SegValue< 0.9 & SegValue> -1.2)
    {
      
      
      Segment1CopyTest<-NWindows*Median1CopyTest
      
 
      IntervalStep1<-seq(1,160,by=1)/100
      
      
      #### Inferring Copy Number 3 ####
      if (SegValue>=0)
      {
        Copy3LikelihoodVec<-c()
        Copy3LikelihoodVecConf2<-c()
        Copy3LikelihoodVecConf1<-c()
        
        for (ii in 1:length(IntervalStep1))
        {
          
          MuStep1<-2*Segment1CopyTest+((Segment1CopyTest)*IntervalStep1[ii])
          Copy3LikelihoodVecConf2[ii]<-dpois(round(TotalReadsSeg+1.96*sqrt(TotalReadsSeg)), MuStep1, log = FALSE)
          Copy3LikelihoodVecConf1[ii]<-dpois(round(TotalReadsSeg-1.96*sqrt(TotalReadsSeg)), MuStep1, log = FALSE)
          
          Copy3LikelihoodVec[ii]<-dpois(round(TotalReadsSeg), MuStep1, log = FALSE)
          
        }
        indLikeliMax<-which.max(Copy3LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy3LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy3LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Duplication",3,FractionEstimate,FractionConf,SegValue))
      }
      
      #### Inferring Copy Number 1 ####
      if (SegValue<0)
      {
        Copy1LikelihoodVec<-c()
        Copy1LikelihoodVecConf2<-c()
        Copy1LikelihoodVecConf1<-c()
        for (ii in 1:length(IntervalStep1))
        {
          
          MuStep1<-2*Segment1CopyTest-((Segment1CopyTest)*IntervalStep1[ii])
          Copy1LikelihoodVecConf2[ii]<-dpois(round(TotalReadsSeg+1.96*sqrt(TotalReadsSeg)), MuStep1, log = FALSE)
          Copy1LikelihoodVecConf1[ii]<-dpois(round(TotalReadsSeg-1.96*sqrt(TotalReadsSeg)), MuStep1, log = FALSE)
          
          Copy1LikelihoodVec[ii]<-dpois(round(TotalReadsSeg), MuStep1, log = FALSE)
        }
        indLikeliMax<-which.max(Copy1LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy1LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy1LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Deletion",1,FractionEstimate,FractionConf,SegValue))
        
      }
      
    }
    
  }
  MatrixOut
}


FractionCallTwoSample<-function(TotalPredBreak,DataSeqTestChrom,DataSeqRefChrom,DataSeq,Median1CopyTest,Median1CopyRef,Pos,chr)
{
  library(skellam)
  MatrixOut<-c()
  for (ll in 1:(length(TotalPredBreak)-1))
  {
    indSeg<-c((TotalPredBreak[ll]+1):TotalPredBreak[ll+1])
    NWindows<-length(indSeg)
    SegValue<-median(DataSeq[indSeg])
 
    if (SegValue>= 0.9)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Multi-Copy Duplication",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    if (SegValue<= -1.2)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Double Deletion",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    
    if (SegValue< 0.9 & SegValue> -1.2)
    {
      
      Segment1CopyTest<-NWindows*Median1CopyTest
      Segment1CopyRef<-NWindows*Median1CopyRef
      
      DiffRC<-round(sum(DataSeqTestChrom[indSeg])-sum(DataSeqRefChrom[indSeg]))
      SumRC<-round(sum(DataSeqTestChrom[indSeg])+sum(DataSeqRefChrom[indSeg]))
      
      
      
      IntervalStep1<-seq(1,160,by=1)/100
      
      
      #### Inferring Copy Number 3 ####
      if (SegValue>=0)
      {
        Copy3LikelihoodVec<-c()
        Copy3LikelihoodVecConf2<-c()
        Copy3LikelihoodVecConf1<-c()
        
        for (ii in 1:length(IntervalStep1))
        {
          
          MuStep1<-2*Segment1CopyTest+((Segment1CopyTest)*IntervalStep1[ii])
          Copy3LikelihoodVecConf2[ii]<-dskellam.sp(DiffRC+1.96*sqrt(SumRC), MuStep1, 2*Segment1CopyRef, log = FALSE)
          Copy3LikelihoodVecConf1[ii]<-dskellam.sp(DiffRC-1.96*sqrt(SumRC), MuStep1, 2*Segment1CopyRef, log = FALSE)
          
          Copy3LikelihoodVec[ii]<-dskellam.sp(DiffRC, MuStep1, 2*Segment1CopyRef, log = FALSE)
          
        }
        indLikeliMax<-which.max(Copy3LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy3LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy3LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Duplication",3,FractionEstimate,FractionConf,SegValue))
      }
      
      #### Inferring Copy Number 1 ####
      if (SegValue<0)
      {
        Copy1LikelihoodVec<-c()
        Copy1LikelihoodVecConf2<-c()
        Copy1LikelihoodVecConf1<-c()
        for (ii in 1:length(IntervalStep1))
        {
          
          MuStep1<-2*Segment1CopyTest-((Segment1CopyTest)*IntervalStep1[ii])
          Copy1LikelihoodVecConf2[ii]<-dskellam.sp(DiffRC+1.96*sqrt(SumRC), MuStep1, 2*Segment1CopyRef, log = FALSE)
          Copy1LikelihoodVecConf1[ii]<-dskellam.sp(DiffRC-1.96*sqrt(SumRC), MuStep1, 2*Segment1CopyRef, log = FALSE)
          
          Copy1LikelihoodVec[ii]<-dskellam.sp(DiffRC, MuStep1,2*Segment1CopyRef, log = FALSE)
        }
        indLikeliMax<-which.max(Copy1LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy1LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy1LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Deletion",1,FractionEstimate,FractionConf,SegValue))
        
      }
      
    }
    
  }
  MatrixOut
}




###### Funzione per Salvare i risultati di FastCall per ogni Regione del target in formato VCF 4.0 ####
VCFRegionCreate<-function(Assembly,DataFolder,ExpLabelOut,TargetFolder,SummaryData,MetaData,out)
{
  StringFormat<-'##fileformat=VCFv4.0'
  StringDate<-paste("##fileDate=",format(Sys.time(), "%Y%d%m"),sep="")
  StringSource<-'##source=Xcavator0.x'
  StringReference<-paste("##reference=",Assembly,sep="")
  StringAssembly<-'##assembly=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta'
  StringINFO1<-'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
  StringINFO2<-'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  StringINFO3<-'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
  StringINFO4<-'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
  StringFORMAT1<-'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  StringFORMAT2<-'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
  StringFORMAT3<-'##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Copy number genotype fraction for imprecise events">'
  StringFORMAT4<-'##FORMAT=<ID=FCL,Number=1,Type=Float,Description="Label Inferred by FastCall algorithm">'
  StringFORMAT5<-'##FORMAT=<ID=FCP,Number=1,Type=Float,Description="Posterior Probability inferred by FastCall algorithm">'
  StringALT1<-'##ALT=<ID=CNV,Description="Copy number variable region">'
  HeaderMat<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
  Header<-c(StringFormat,StringDate,StringSource,StringReference,StringAssembly,StringINFO1,StringINFO2,StringINFO3,StringINFO4,StringALT1,StringFORMAT1,StringFORMAT2,StringFORMAT3,StringFORMAT4,StringFORMAT5)
  FormatField<-"GT:CN:CNF:FCL:FCP"
  HeadMatSample<-c(HeaderMat,ExpLabelOut)
  
  
  indSig<-which(out[,1]!=0)
  
  
  
  
  if (length(indSig)!=0)
  {
    
    
    outSig<-out[indSig,]
    P0Sig<-P0[indSig,]
    SummaryDataSig<-SummaryData[indSig,]
    CNFSig<-round(2*(2^SummaryDataSig[,4]),digit=2)
    CNSig<-round(CNFSig)
    
    
    MatSig<-matrix(0,nrow=length(indSig),ncol=7)
    
    
    for (j in 1:length(indSig))
    {
      
      indS<-SummaryDataSig[j,2]
      indE<-SummaryDataSig[j,3]
      MatSig[j,]<-cbind(MetaData[indS,1],MetaData[indS,3],MetaData[indE,4],CNFSig[j],CNSig[j],outSig[j,1],round(outSig[j,2],digit=2))
    }
    
    
    ChromoSig<-MatSig[,1]
    StartSig<-MatSig[,2]
    EndSig<-MatSig[,3]
    CNFSig<-MatSig[,4]
    CNSig<-MatSig[,5]
    CallSig<-MatSig[,6]
    ProbSig<-MatSig[,7]
    
    
    
    indEnde<-grep("e",EndSig)
    EndSige<-EndSig[indEnde]
    indStarte<-grep("e",StartSig)
    StartSige<-StartSig[indStarte]
    
    options(scipen = 10)
    
    StartSig[indStarte]<-as.character(as.numeric(StartSige))
    EndSig[indEnde]<-as.character(as.numeric(EndSige))
    L<-as.character(as.numeric(EndSig)-as.numeric(StartSig)+1)
    
    MatOut<-matrix(NA,ncol=10,nrow=length(ChromoSig))
    ChromoSigU<-unique(ChromoSig)
    
    indStart<-1
    for (i in 1:length(ChromoSigU))
    {
      
      FRBIn<-file.path(TargetFolder,"FRB",paste("FRB.",ChromoSigU[i],".RData",sep=""))
      load(FRBIn)
      FRBPos<-as.character(FRBData[,1])
      FRBRef<-as.character(FRBData[,2])
      
      indC<-which(ChromoSig==ChromoSigU[i])
      ChromoSigC<-ChromoSig[indC]
      EndSigC<-EndSig[indC]
      StartSigC<-StartSig[indC]
      CNSigC<-CNSig[indC]
      CNFSigC<-CNFSig[indC]
      CallSigC<-CallSig[indC]
      ProbSigC<-ProbSig[indC]
      LC<-L[indC]
      RefC<-FRBRef[which(!is.na(match(FRBPos,StartSigC)))]
      AltC<-rep("<CNV>",length(RefC))
      indEnd<-indStart+length(indC)-1
      InfoField<-paste("IMPRECISE;SVTYPE=CNV;END=",EndSigC,";SVLEN=",LC,";",sep="")
      GenoField<-paste("1/1:",CNSigC,":",CNFSigC,":",CallSigC,":",ProbSigC,sep="")
      MatOut[c(indStart:indEnd),]<-cbind(ChromoSigC,StartSigC,".",RefC,AltC,".","PASS",InfoField,FormatField,GenoField)
      indStart<-indEnd+1
    }
  }
  if (length(indSig)==0)
  {
    MatOut<-c()
  }
  MatOut<-rbind(HeadMatSample,MatOut)
  
  FileOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("XcavatorRegionCall_",ExpLabelOut,".vcf",sep=""))
  
  zz <- file(FileOut, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  
}








