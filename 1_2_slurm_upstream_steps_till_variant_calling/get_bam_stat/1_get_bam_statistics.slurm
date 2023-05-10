#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --export="NONE"
#SBATCH --job-name="1_get_bam_statistics.pbs"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="2"
#SBATCH --ntasks="2"
#SBATCH --output="xx/%x.o%A"
#SBATCH --time="10:00:00"

# sbatch --cluster=genius 1_get_bam_statistics.slurm

module load SAMtools/1.9-GCC-6.4.0-2.28

SAMPLE="xx"
bam_path="xx"
cd $bam_path

# samtools flagstat
samtools_flagstat=`samtools flagstat $SAMPLE.trimmed.sorted.bam`

# get num of chimeric reads
num_read_chimera=`samtools view -f 0x800 $SAMPLE.trimmed.sorted.bam|cut -f1|sort|uniq|wc -l`

# get num of multimapped reads
num_read_multimapping=`samtools view -f 256 $SAMPLE.trimmed.sorted.bam|cut -f1|sort|uniq|wc -l`

# get num of unmapped reads （flag为4）
num_read_unmapped=`samtools view -f 4 $SAMPLE.trimmed.sorted.bam|wc -l`

echo -e "samtools_flagstat\n" $samtools_flagstat "\n\nnum_read_chimera" $num_read_chimera "\nnum_read_multimapping" $num_read_multimapping "\nnum_read_unmapped" $num_read_unmapped > $SAMPLE.bam.stats.txt
