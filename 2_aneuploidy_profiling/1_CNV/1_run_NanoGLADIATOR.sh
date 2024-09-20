# script1 ran already!!! start from script 2!!!!!!!!
# required software and packages


module use /apps/leuven/rocky8/skylake/2018a/modules/all
module load R/4.0.2-foss-2018a-bare 
module load SAMtools/1.9-GCC-6.4.0-2.28
# module load Python/2.7.14-foss-2018a (use python3 instead)
module load BEDTools/2.27.1-intel-2018a

cd /vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/NanoGLADIATOR_1.0

# revise NanoGLADIATORPrepare.pl  BamFilePrepare.txt BamFileAnalysis.txt
# run script 2,3

# -----script2: NanoGLADIATORPrepare.pl--------

# change file 1: path to an input text file: BamFilePrepare.txt
# bampath outpath sample_name

 nohup perl NanoGLADIATORPrepare.pl BamFilePrepare.txt --processors 6 --target hg38w1000000 --assembly hg38 & 

# -------script3: NanoGLADIATORAnalysis.pl-------
# change file 2: BamFileAnalysis.txt (last 2 columns same as ast 2 columns in above BamFilePrepare.txt)

# !! revise --output in command below
nohup perl NanoGLADIATORAnalysis.pl BamFileAnalysis.txt --processors 6 --target hg38w1000000 --assembly hg38 --output /lustre1/project/stg_00019/research/yan/nanopore_data/01_03_human_PGT_families/CNV/hg38w1000000/ONT1-E02 --mode nocontrol &

