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


###########
COMMANDS START
#############

bowtie2 -x /home/ubuntu/ambi/genome/PA7Genome -1 qualityreads/R_output_R1_forward_paired.fq.gz -2 qualityreads/R_output_R2_reverse_paired.fq.gz -S alignment/R_reads_aligned.sam\
&& samtools view -b -S -o alignment/R_reads_aligned.bam alignment/R_reads_aligned.sam\
&&samtools sort -o alignment/R_reads_aligned_sorted.bam alignment/R_reads_aligned.bam\
&& samtools index alignment/R_reads_aligned_sorted.bam\
&& samtools mpileup -g -f /home/ubuntu/ambi/genome/PA7Genome.fna alignment/R_reads_aligned_sorted.bam > variants/R_variants.bcf\
&& bcftools call -c -v variants/R_variants.bcf > variants/R_variants.vcf

freebayes -f /home/ubuntu/ambi/genome/normalized/PA7Genome3.fna alignment/R_reads_aligned_sorted.bam > variants/R_Freebayes_variants.vcf

##########
START WORKFLOW 1
##########

prefix="R"
path1=/home/ubuntu/ambi

bowtie2 -x $path1/genome/denovo/WT_scaffolds_denovo_filtered -1 qualityreads/"$prefix"_output_R1_forward_paired.fq.gz -2 qualityreads/"$prefix"_output_R2_reverse_paired.fq.gz -S denovo/alignment/"$prefix"_reads_aligned_denovo.sam\
&&samtools view -b -S -o denovo/alignment/"$prefix"_reads_aligned_denovo.bam denovo/alignment/"$prefix"_reads_aligned_denovo.sam\
&&samtools sort -o denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.bam denovo/alignment/"$prefix"_reads_aligned_denovo.bam\
&&samtools index denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.bam\
&& freebayes -f $path1/genome/denovo/WT_scaffolds_denovo_filtered.fasta denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.bam > variants/"$prefix"_variants_freebayes_denovo.vcf\

##########
END WORKFLOW 1
##########

##########
START WORKFLOW 2
##########

prefix="R"
path1=/home/ubuntu/ambi
pathpicard=/home/ubuntu/ambi/software/Picard/

java -Xms4g -jar $pathpicard/picard.jar MarkDuplicates INPUT=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.bam OUTPUT=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam METRICS_FILE=denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.txt2 REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT\
&&samtools index denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam\
&& freebayes -f $path1/genome/denovo/WT_scaffolds_denovo_filtered.fasta denovo/alignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam > variants/"$prefix"_variants_freebayes_denovo_rmdup.vcf\



##########
END WORKFLOW 2
##########


END


##########
#START OF BACKUP OF MAIN WORKFLOW 3 APRIL 2017
##########

prefix="R"
#path1=/home/ubuntu/ambi
#path softwares
pathPicard="/home/ubuntu/ambi/software/Picard/"
pathFastqc="/home/ubuntu/ambi/software/fastqc_v0.11.5/FastQC/"
pathTrimmomaticAdapters="/home/ubuntu/ambi/software/Trimmomatic-0.36/adapters"
pathVarscan="/home/ubuntu/ambi/software/varscan/"


#path raw fastq files containing reads
#path2=qualitytrimreads
pathQualityReads="qualitytrimreads"

#readsFile1="B_S60_L001_R1_001_val_1.fq.gz"
#readsFile2="B_S60_L001_R2_001_val_2.fq.gz"

readsFile1="/home/ubuntu/ambi/scvB_Rev/R_S62_L001_R1_001_val_1.fq.gz"
readsFile2="/home/ubuntu/ambi/scvB_Rev/R_S62_L001_R2_001_val_2.fq.gz"


pathAlignment="alignment"

#path genome
pathGenome=/home/ubuntu/ambi/genome/denovo
pathVariants=variants

mkdir $pathQualityReads
mkdir "$pathQualityReads"/fastqcoutput_"$prefix"
mkdir $pathAlignment
mkdir $pathVariants

pathQualityReadsForwardPaired="$pathQualityReads"/"$prefix"_output_forward_paired.fq.gz
pathQualityReadsForwardUnPaired="$pathQualityReads"/"$prefix"_output_forward_unpaired.fq.gz

pathQualityReadsReversePaired="$pathQualityReads"/"$prefix"_output_reverse_paired.fq.gz
pathQualityReadsReverseUnPaired="$pathQualityReads"/"$prefix"_output_reverse_unpaired.fq.gz


#trim reads
trimmomatic PE -phred33 "$readsFile1" "$readsFile2" $pathQualityReadsForwardPaired $pathQualityReadsForwardUnPaired $pathQualityReadsReversePaired $pathQualityReadsReverseUnPaired  ILLUMINACLIP:"$pathTrimmomaticAdapters"/TruSeq3-PE.fa:2:30:10 CROP:160 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36\
&&echo "generate fastqc report"\
&&"$pathFastqc"/fastqc $pathQualityReadsForwardPaired $pathQualityReadsReversePaired --outdir="$pathQualityReads"/fastqcoutput_"$prefix"/\
&&echo "map using bowtie2"\
&&bowtie2 -x $pathGenome/WT_scaffolds_denovo_filtered -1 $pathQualityReadsForwardPaired -2 $pathQualityReadsReversePaired -S $pathAlignment/"$prefix"_reads_aligned_denovo.sam\
&&echo "convert sam to bam using samtools view"\
&&samtools view -b -S -o $pathAlignment/"$prefix"_reads_aligned_denovo.bam $pathAlignment/"$prefix"_reads_aligned_denovo.sam\
&&echo "sort bam using samtools sort"\
&&samtools sort -o $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.bam $pathAlignment/"$prefix"_reads_aligned_denovo.bam\
&&echo "remove duplicates"\
&&java -Xms4g -jar $pathPicard/picard.jar MarkDuplicates INPUT= $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.bam OUTPUT= $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam METRICS_FILE= $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.txt2 REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT\
&&echo "generate index of bam file to view it in igv"\
&&samtools index $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam\
&&echo "generate fastqc report post remove duplicate"\
&&"$pathFastqc"/fastqc $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam --outdir="$pathQualityReads"/fastqcoutput_"$prefix"/\
&&echo "generate vcf using freebayes"\
&& freebayes -f $pathGenome/WT_scaffolds_denovo_filtered.fasta $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam > $pathVariants/"$prefix"_variants_freebayes_denovo_rmdup.vcf\
&&echo "generate mpileup"\
&&samtools mpileup -f $pathGenome/WT_scaffolds_denovo_filtered.fasta $pathAlignment/"$prefix"_reads_aligned_denovo_sorted.rmdup.bam > "$prefix".mpileup\
&&echo "generate vcf using varscan and mpileup using varscan mpileup2cns"\
&&java -jar $pathVarscan/VarScan.v2.3.9.jar mpileup2cns "$prefix".mpileup --min-coverage 40 --min-reads2 20 --output-vcf 1 --variants > "$prefix".vcf


##########
#END OF BACKUP OF MAIN WORKFLOW 3 APRIL 2017
##########

