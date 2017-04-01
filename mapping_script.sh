#!/bin/bash
: <<'END'

Tutorials
1. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs
2. http://biobits.org/samtools_primer.html

Indexing a reference genome
bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus

Aligning reads Paired-end example
$BT2_HOME/bowtie2 -x lambda_virus -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam

Use samtools view to convert the SAM file into a BAM file. BAM is a the binary format corresponding to the SAM text format. Run:

samtools view -bS eg2.sam > eg2.bam
Use samtools sort to convert the BAM file to a sorted BAM file.

samtools sort eg2.bam eg2.sorted

samtools mpileup -g -f genomes/NC_008253.fna alignments/sim_reads_aligned.sorted.bam > variants/sim_variants.bcf

W_S63_L001_R1_001_val_1.fq.gz

W_S63_L001_R2_001_val_2.fq.gz


trimmomatic PE -phred33 W_S63_L001_R1_001_val_1.fq.gz W_S63_L001_R2_001_val_2.fq.gz WT_output_forward_paired.fq.gz WT_output_forward_unpaired.fq.gz WT_output_reverse_paired.fq.gz WT_output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic SE sra_data.fastq.gz sra_data.fastq_trim.gz SLIDINGWINDOW:4:20 MINLEN:20


#############
Set of commands

$BT2_HOME/bowtie2 -x lambda_virus -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam

samtools view -b -S -o alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sam
samtools sort alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sorted
samtools index alignments/sim_reads_aligned.sorted.bam
samtools mpileup -g -f genomes/NC_008253.fna alignments/sim_reads_aligned.sorted.bam > variants/sim_variants.bcf
bcftools call -c -v variants/sim_variants.bcf > variants/sim_variants.vcf

samtools view -bS eg2.sam > eg2.bam
samtools sort eg2.bam eg2.sorted

samtools mpileup -g -f genomes/NC_008253.fna alignments/sim_reads_aligned.sorted.bam > variants/sim_variants.bcf

END

bowtie2 -x ../genome/PA7Genome -1 qualityreads/WT_output_forward_paired.fq.gz -2 qualityreads/WT_output_reverse_paired.fq.gz -S alignment/WT_reads_aligned.sam
samtools view -b -S -o alignment/WT_reads_aligned.bam alignment/WT_reads_aligned.sam
samtools sort alignment/WT_reads_aligned.bam -o alignment/WT_reads_aligned_sorted.bam
samtools index alignment/WT_reads_aligned_sorted.bam
samtools mpileup -g -f ../genome/PA7Genome.fna alignment/WT_reads_aligned_sorted.bam > variants/WT_variants.bcf
bcftools call -c -v variants/WT_variants.bcf > variants/WT_variants.vcf
