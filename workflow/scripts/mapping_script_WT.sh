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

############
COMMANDS SRARTS
#############

bowtie2 -x /home/ubuntu/ambi/genome/PA7Genome -1 qualityreads/WT_output_forward_paired.fq.gz -2 qualityreads/WT_output_reverse_paired.fq.gz -S alignment/WT_reads_aligned.sam\
&&samtools view -b -S -o alignment/WT_reads_aligned.bam alignment/WT_reads_aligned.sam\
&&samtools sort -o alignment/WT_reads_aligned_sorted.bam alignment/WT_reads_aligned.bam\
&&samtools index alignment/WT_reads_aligned_sorted.bam\
&&samtools mpileup -g -f /home/ubuntu/ambi/genome/PA7Genome.fna alignment/WT_reads_aligned_sorted.bam > variants/WT_variants.bcf\
&&bcftools call -c -v variants/WT_variants.bcf > variants/WT_variants.vcf


##########
START WORKFLOW 1
##########

bowtie2 -x /home/ubuntu/ambi/genome/denovo/WT_scaffolds_denovo_filtered -1 qualityreads/WT_output_forward_paired.fq.gz -2 qualityreads/WT_output_reverse_paired.fq.gz -S denovo/alignment/WT_reads_aligned_denovo.sam\
&&samtools view -b -S -o denovo/alignment/WT_reads_aligned_denovo.bam denovo/alignment/WT_reads_aligned_denovo.sam\
&&samtools sort -o denovo/alignment/WT_reads_aligned_denovo_sorted.bam denovo/alignment/WT_reads_aligned_denovo.bam\
&&samtools index denovo/alignment/WT_reads_aligned_denovo_sorted.bam\
&& freebayes -f /home/ubuntu/ambi/genome/denovo/WT_scaffolds_denovo_filtered.fasta denovo/alignment/WT_reads_aligned_denovo_sorted.bam > variants/WT_variants_freebayes_denovo.vcf\
&&samtools mpileup -g -f /home/ubuntu/ambi/genome/denovo/WT_scaffolds_denovo_filtered.fasta denovo/alignment/WT_reads_aligned_denovo_sorted.bam > variants/WT_variants_denovo.bcf\
&&bcftools call -c -v variants/WT_variants_denovo.bcf > variants/WT_variants_denovo.vcf

##########
END WORKFLOW 1
##########

END

prefix="WT"
path1=/home/ubuntu/ambi
pathpicard=/home/ubuntu/ambi/software/Picard/

java -Xms4g -jar $pathpicard/picard.jar MarkDuplicates INPUT=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.bam OUTPUT=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam METRICS_FILE=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.txt2 REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT\
&&samtools index denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam\
&& freebayes -f $path1/genome/denovo/WT_scaffolds_denovo_filtered.fasta denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam > variants/"$prefix"_variants_freebayes_denovo_rmdup.vcf\


