#!/bin/bash
#$ -N align
#$ -pe thread 16
#$ -cwd
#$ -j y
#$ -M soumya.negi@biogen.com
#$ -m e
#$ -o align.log

source /etc/bashrc
module load bwa/0.7.17
module load samtools/1.10

reference_V79="/camhpc/ngs/projects/TST11724/AAV_genome/V79/V79_DN_145_pERG.00249_SOD1_miRE_7_lambdaKanlambda.fa"

# Building reference
#bwa index $reference_V79/V79_DN_145_pERG.00249_SOD1_miRE_7_lambdaKanlambda.fa

fastq_dir="/camhpc/ngs/projects/TST11724/DNANexus/unmapped_fastq"
samples="/camhpc/ngs/projects/TST11724/scripts/sample_sheet.txt"
bam_dir="/camhpc/ngs/projects/TST11724/alignment_mRNAseq/From_Unmapped_reads/complete_bam"
cd $bam_dir
date
while read line; do
	echo $line
	file1=$fastq_dir/TST11580_"$line"_unmapped_R1.fastq
	file2=$fastq_dir/TST11580_"$line"_unmapped_R2.fastq
	echo $file1
	echo $file2
	
	bam_file=$bam_dir/"$line".unmapped.mRNA.V79_AAV-amiR-SOD1.bam
	echo "Running bwa on unmapped $line"
	bwa mem -P \
	-t 8 \
	$reference_V79 \
	$file1 $file2 | samtools view -bS -@8 | samtools sort -@8 -o $bam_file
	
	echo "Indexing $line"
	samtools index -@8 $bam_file 
	
	echo "Sorting $line"
	samtools sort -@8 $bam_file 
	
	echo "Determining coverage of $line"
	samtools coverage -A -o $bam_dir/"$line".unmapped.mRNA.AAV.coverage_histogram.txt $bam_file 
	samtools coverage -o $bam_dir/"$line".unmapped.mRNA.AAV.coverage.txt $bam_file 
	
	echo "Running flagstat on $line"
	samtools flagstat -@ 8 $bam_file  > $bam_dir/"$line".unmapped.mRNA.AAV.mapping.stat.txt

done < $samples

# Code below splits sense and antisense transcripts into separate bams
sense_antisense_dir="/camhpc/ngs/projects/TST11724/alignment_mRNAseq/From_Unmapped_reads/sense_antisense"
cd $sense_antisense_dir
while read line; do
	echo $line
	file=$(ls $bam_dir/$line*.bam)
	echo $file
	
	# Forward strand:
	samtools view -@8 -b -f 128 -F 16 $file > fwd1.bam
	samtools index -@8 fwd1.bam

	samtools view -@8 -b -f 80 $file > fwd2.bam
	samtools index -@8 fwd2.bam
	
	sense_bam=$sense_antisense_dir/"$line".unmapped.mRNA.V79_AAV-amiR-SOD1_SENSE.bam
	samtools merge -f $sense_bam fwd1.bam fwd2.bam
	samtools index -@8 $sense_bam
	echo "Running indexing on sense $line"
	samtools sort -@8 $sense_bam
	echo "Running coverage on sense $line"
	samtools coverage -A -o $sense_bam.coverage_histogram.txt $sense_bam 
	samtools coverage -o $sense_bam.coverage.txt $sense_bam 
	echo "Running flagstat on sense $line"
	samtools flagstat -@ 8 $sense_bam  > $sense_bam.mapping.stat.txt
	
	# Reverse strand

	samtools view -@8 -b -f 144 $file > rev1.bam
	samtools index -@8 rev1.bam

	samtools view -@8 -b -f 64 -F 16 $file > rev2.bam
	samtools index -@8 rev2.bam

	antisense_bam=$sense_antisense_dir/"$line".unmapped.mRNA.V79_AAV-amiR-SOD1_ANTISENSE.bam
	samtools merge -@8 -f $antisense_bam rev1.bam rev2.bam
	samtools index -@8 $antisense_bam
	echo "Running indexing on antisense $line"
	samtools sort -@8 $antisense_bam
	echo "Running coverage on antisense $line"
	samtools coverage -A -o $antisense_bam.coverage_histogram.txt $antisense_bam 
	samtools coverage -o $antisense_bam.coverage.txt $antisense_bam 
	echo "Running flagstat on antisense $line"
	samtools flagstat -@ 8 $antisense_bam  > $antisense_bam.mapping.stat.txt
	
	rm fwd1.bam fwd1.bam.bai fwd2.bam fwd2.bam.bai rev1.bam rev1.bam.bai rev2.bam rev2.bam.bai

done < $samples
