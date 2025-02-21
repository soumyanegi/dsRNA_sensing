#!/bin/bash
#$ -N extract_unaligned_reads
#$ -pe thread 16
#$ -cwd
#$ -j y
#$ -M soumya.negi@biogen.com
#$ -m e
#$ -o extract_unaligned_reads.log

source /etc/bashrc
module load samtools/1.10
module load bedtools/2.29.2

samples="/camhpc/ngs/projects/TST11724/scripts/sample_sheet.txt"
bam_dir="/camhpc/ngs/projects/TST11724/DNANexus"
fastq_output="/camhpc/ngs/projects/TST11724/DNANexus/unmapped_fastq"

cd $fastq_output

date
while read line; do
	echo $line
	file=$bam_dir/TST11580_"$line".genome.sorted.bam
	echo $file
	
	echo "Extracting unmapped reads on $line"
	
	samtools view -@ 8 -u  -f 4 -F264 $file  > "$line".tmps1.bam
	samtools view -@ 8 -u -f 8 -F 260 $file  > "$line".tmps2.bam
	samtools view -@ 8 -u -f 12 -F 256 $file > "$line".tmps3.bam
	
	samtools merge -@ 8 -u - "$line".tmps1.bam "$line".tmps2.bam "$line".tmps3.bam | samtools sort -@ 8 -n -o TST11580_"$line"_unmapped.bam
	
	#samtools view -@ 8 -S -b TST11580_"$line"_unmapped.sam > TST11580_"$line"_unmapped.bam
	bamToFastq -i TST11580_"$line"_unmapped.bam -fq TST11580_"$line"_unmapped_R1.fastq -fq2 TST11580_"$line"_unmapped_R2.fastq

 	rm "$line".tmps1.bam "$line".tmps2.bam "$line".tmps3.bam
done < $samples