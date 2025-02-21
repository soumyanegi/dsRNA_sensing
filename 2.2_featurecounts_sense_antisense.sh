#!/bin/bash
#$ -N featurecount
#$ -pe thread 16
#$ -cwd
#$ -j y
#$ -M soumya.negi@biogen.com
#$ -m e
#$ -o featurecount.log


# Code adapted from Fergals: /camhpc/ngs/projects/TST11147/look_at_featurecounts/featurecounts_v27_v31only.sh
source /etc/bashrc
module load subread/1.6.4

samples="/camhpc/ngs/projects/TST11724/scripts/sample_sheet_pilot.txt"
bam_dir="/camhpc/ngs/projects/TST11724/alignment_mRNAseq/DNANexus_20210226215854/sense_antisense"
bam_dir2="/camhpc/ngs/projects/TST11724/DNANexus/20210317213939/sort_bams"
GTF_FILE="/camhpc/ngs/projects/TST11724/MacFas_AAV_merged/AAV_sense_antisense_directionality/macaca_fascicularis_5.0_V79.AAV.amiR.SOD1_w_directionality.gtf"
output="/camhpc/ngs/projects/TST11724/counts_mRNAseq/pairs"

cd $bam_dir

date
while read line; do
	
	
	STRAND=2
	
	# --minOverlap 1 is absolutely essential to assign counts to the smaller ITR fragments. DNANexus uses --minOverlap 25, which doesn't work
	FEATURECOUNTS_OVERLAP="--minOverlap 1"
	PAIRS="-p -B"
	
	bam=$bam_dir2/TST11580_"$line".genome.sorted.bam
	featureCounts -O -T 16 -F GTF -a $GTF_FILE -t exon -g gene_id -s $STRAND $PAIRS $FEATURECOUNTS_OVERLAP -o $output/"$line".mRNA.V79_AAV-amiR-SOD1.s1.ALL.counts.txt $bam
	
	bamqis_1=$bam_dir/"$line".mRNA.V79_AAV-amiR-SOD1_SENSE.bam
	#Section 6.2.6
	#http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf
	# By default, featureCounts does not count multi-overlapping reads. 
	# Users can specify the -O option (set allowMultiOverlap to TRUE in R) to fully count them for each overlapping metafeature/feature 
	# (each overlapping meta-feature/feature receives a count of 1 from a read), or specify both -O and -fraction
	# options (set both allowMultiOverlap and fraction to TRUE in R) to assign a fractional count to each overlapping 
	# meta-feature/feature (each overlapping meta-feature/feature receives a count of 1/y from a read where y is the total number of meta-features/features overlapping with the read).
	featureCounts -O -T 16 -F GTF -a $GTF_FILE -t exon -g gene_id -s $STRAND $PAIRS $FEATURECOUNTS_OVERLAP -o $output/"$line".mRNA.V79_AAV-amiR-SOD1.s1.SENSE.counts.txt $bamqis_1
	
	bamqis_2=$bam_dir/"$line".mRNA.V79_AAV-amiR-SOD1_ANTISENSE.bam
	
	featureCounts -O -T 16 -F GTF -a $GTF_FILE -t exon -g gene_id -s $STRAND $PAIRS $FEATURECOUNTS_OVERLAP -o $output/"$line".mRNA.V79_AAV-amiR-SOD1.s1.ANTISENSE.counts.txt $bamqis_2

done < $samples
