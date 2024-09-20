#!/bin/bash

# chmod 755 4_nanoplot_fastq_trimmed_get_stat.sh
# ./4.1_nanoplot_fastq_trimmed_get_stat.sh


Sample=""
wkdir=""
dir=$wkdir/4_nanoplot_fastq_trimmed

cd $dir

Mean_read_length=`cat NanoStats.txt|grep "Mean read length"|awk '{print $NF}'`
Mean_read_quality=`cat NanoStats.txt|grep "Mean read quality"|awk '{print $NF}'`
Median_read_length=`cat NanoStats.txt|grep "Median read length"|awk '{print $NF}'`
Median_read_quality=`cat NanoStats.txt|grep "Median read quality"|awk '{print $NF}'`
Number_of_reads=`cat NanoStats.txt|grep "Number of reads"|awk '{print $NF}'`
Read_length_N50=`cat NanoStats.txt|grep "Read length N50"|awk '{print $NF}'`
Stdev_read_length=`cat NanoStats.txt|grep "STDEV read length"|awk '{print $NF}'`	
Total_bases=`cat NanoStats.txt|grep "Total bases"|awk '{print $NF}'`
Q5=`cat NanoStats.txt|grep "Q5"|awk -F[\(\)] '{print $2}'`
Q7=`cat NanoStats.txt|grep "Q7"|awk -F[\(\)] '{print $2}'`
Q10=`cat NanoStats.txt|grep "Q10"|awk -F[\(\)] '{print $2}'`	
Q12=`cat NanoStats.txt|grep "Q12"|awk -F[\(\)] '{print $2}'`	
Q15=`cat NanoStats.txt|grep "Q15"|awk -F[\(\)] '{print $2}'`

echo -e Sample"\t"Fastq_file"\t"Mean_read_length"\t"Mean_read_quality"\t"Median_read_length"\t"Median_read_quality"\t"Number_of_reads"\t"Read_length_N50"\t"Stdev_read_length"\t"Total_bases"\t"Q5"\t"Q7"\t"Q10"\t"Q12"\t"Q15>clean_stat.txt
echo -e $Sample"\t"trimmed"\t"$Mean_read_length"\t"$Mean_read_quality"\t"$Median_read_length"\t"$Median_read_quality"\t"$Number_of_reads"\t"$Read_length_N50"\t"$Stdev_read_length"\t"$Total_bases"\t"$Q5"\t"$Q7"\t"$Q10"\t"$Q12"\t"$Q15>>clean_stat.txt