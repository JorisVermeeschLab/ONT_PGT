#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --get-user-env="60L"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="2"
#SBATCH --ntasks="2"
#SBATCH --time="20:00:00"

# sbatch --cluster=genius 5.3_count_depthofcoverage.slurm

SAMPLE=""
file_path="xx/5_mapped"
cd $file_path
file=$SAMPLE.samtools_depth.txt

# calculate average depth of coverage
total_bp=`cat $file|wc -l`
sum_of_depth_all_bp=`cat $file|awk '{sum+=$3}END{print sum}'`
average_depth=`echo -e $sum_of_depth_all_bp"\t"$total_bp|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`
num_bp_depth_0=`cat $file|awk '{if($3==0)print}'|wc -l`
perc_genome_covered=`echo -e $num_bp_depth_0"\t"$total_bp|awk 'BEGIN{OFMT = "%.4f"}{print (1-$1/$2)}'`

echo -e average_depth"\t"$average_depth"\n"perc_genome_covered"\t"$perc_genome_covered>>depth_summay.txt
