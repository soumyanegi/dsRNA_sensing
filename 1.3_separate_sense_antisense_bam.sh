#!/bin/bash
#$ -N separate_sense_antisense
#$ -pe thread 16
#$ -cwd
#$ -j y
#$ -M soumya.negi@biogen.com
#$ -m e
#$ -o separate_sense_antisense.log

# Code below splits sense and antisense transcripts into separate bams
source /etc/bashrc
module load samtools/1.10

samples="/camhpc/ngs/projects/TST11724/scripts/sample_sheet_pilot.txt"
bam_dir="/camhpc/ngs/projects/TST11724/DNANexus/20210226215854/star_rsem_sample/sort_bams"
sense_antisense_dir="/camhpc/ngs/projects/TST11724/alignment_mRNAseq/DNANexus_20210226215854/sense_antisense"
cd $sense_antisense_dir
while read line; do
	echo $line
	file=$(ls $bam_dir/TST11580_$line*.bam)
	echo $file
	
	# Forward strand:
	samtools view -@8 -b -f 128 -F 16 $file > fwd1.bam
	samtools index -@8 fwd1.bam

	samtools view -@8 -b -f 80 $file > fwd2.bam
	samtools index -@8 fwd2.bam
	
	sense_bam=$sense_antisense_dir/"$line".mRNA.V79_AAV-amiR-SOD1_SENSE.bam
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

	antisense_bam=$sense_antisense_dir/"$line".mRNA.V79_AAV-amiR-SOD1_ANTISENSE.bam
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