#!/bin/bash
#$ -N align_miRNA
#$ -pe thread 16
#$ -cwd
#$ -j y
#$ -M soumya.negi@biogen.com
#$ -m e
#$ -o align_miRNA.log

source /etc/bashrc
module load bwa/0.7.17
module load samtools/1.10

reference_V79="/camhpc/ngs/projects/TST11724/AAV_genome/V79/V79_DN_145_pERG.00249_SOD1_miRE_7_lambdaKanlambda.fa"

# Building reference
#bwa index $reference_V79/V79_DN_145_pERG.00249_SOD1_miRE_7_lambdaKanlambda.fa

fastq_dir="/camhpc/ngs/projects/TST11724/fastq_miRNAseq"
samples="/camhpc/ngs/projects/TST11724/scripts/sample_sheet.txt"
out_dir="/camhpc/ngs/projects/TST11724/alignment_miRNAseq"
date
while read line; do
	echo $line
	file1=$(ls $fastq_dir/$line*L001_R1*)
	file2=$(ls $fastq_dir/$line*L001_R2*)
	echo $file1
	echo $file2
	
	bam_file=$out_dir/"$line".miRNA.V79_AAV-amiR-SOD1.bam
	echo "Running bwa on $line"
	bwa mem -P \
	-t 8 \
	$reference_V79 \
	$file1 $file2 | samtools view -bS -@8 | samtools sort -@8 -o $bam_file
	
	echo "Indexing $line"
	samtools index -@8 $bam_file 
	
	echo "Sorting $line"
	samtools sort -@8 $bam_file 
	
	echo "Determining coverage of $line"
	samtools coverage -A -o $out_dir/"$line".miRNA.AAV.coverage_histogram.txt $bam_file 
	samtools coverage -o $out_dir/"$line".miRNA.AAV.coverage.txt $bam_file 
	
	echo "Running flagstat on $line"
	samtools flagstat -@ 8 $bam_file  > $out_dir/"$line".miRNA.AAV.mapping.stat.txt

done < $samples
