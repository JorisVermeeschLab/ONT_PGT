#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --job-name="3_count_samtools_genome_depthofcoverage.slurm"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="2"
#SBATCH --ntasks="2"
#SBATCH --time="20:00:00"

# sbatch --cluster=genius 3_count_samtools_genome_depthofcoverage.slurm

file_path="xx"
cd $file_path
file="xx.samtools_depth.txt"

# calculate average depth of coverage
total_bp=`cat $file|wc -l`
sum_of_depth_all_bp=`cat $file|awk '{sum+=$3}END{print sum}'`
average_depth=`echo -e $sum_of_depth_all_bp"\t"$total_bp|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`
num_bp_depth_0=`cat $file|awk '{if($3==0)print}'|wc -l`
perc_genome_covered=`echo -e $num_bp_depth_0"\t"$total_bp|awk 'BEGIN{OFMT = "%.4f"}{print (1-$1/$2)}'`

echo -e total_bp"\t"$total_bp"\n"perc_genome_covered"\t"$perc_genome_covered"\n"average_depth"\t"$average_depth>depth_summay.txt

for lim in 5 10 15 20 25 30 100 200
do
num_bp_gt_depth_lim=`cat $file|awk '{if($3>'$lim')print}'|wc -l`
perc_bp_gt_depth_lim=`echo -e $num_bp_gt_depth_lim"\t"$total_bp|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`
echo -e perc_genome_depth_higher_than_$lim"\t"$perc_bp_gt_depth_lim>>depth_summay.txt;
done

for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
	total_bp_chr=`cat $file |awk '{if($1=="'$i'")print}'|wc -l`
	num_bp_depth_0_chr=`cat $file|awk '{if($1=="'$i'")print}'|awk '{if($3==0)print}'|wc -l`
	sum_of_depth_all_bp_chr=`cat $file|awk '{if($1=="'$i'")print}'|awk '{sum+=$3}END{print sum}'`
	perc_chr_covered=`echo -e $num_bp_depth_0_chr"\t"$total_bp_chr|awk 'BEGIN{OFMT = "%.4f"}{print (1-$1/$2)}'`
	average_depth_chr=`echo -e $sum_of_depth_all_bp_chr"\t"$total_bp_chr|awk 'BEGIN{OFMT = "%.4f"}{print $1/$2}'`
	echo -e perc_"$i"_covered"\t"$perc_chr_covered"\n"average_depth_"$i""\t"$average_depth_chr>>depth_summay.txt;
done
