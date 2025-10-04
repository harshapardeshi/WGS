#!/bin/bash

## Before running, make sure these directories exist:
## /mnt/c/Users/Dell/Desktop/yeast/reads
## /mnt/c/Users/Dell/Desktop/yeast/ref
## /mnt/c/Users/Dell/Desktop/yeast/data
## /mnt/c/Users/Dell/Desktop/yeast/aligned_reads
## /mnt/c/Users/Dell/Desktop/yeast/qc_reads
## You can create them using:
## mkdir -p /mnt/c/Users/Dell/Desktop/yeast/{reads,ref,data,aligned_reads,qc_reads}

## set all the directories ##

reads="/mnt/c/Users/Dell/Desktop/yeast/reads"

ref="/mnt/c/Users/Dell/Desktop/yeast/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

data="/mnt/c/Users/Dell/Desktop/yeast/data"      # for all the outputs after the alignement 

aligned_reads="/mnt/c/Users/Dell/Desktop/yeast/aligned_reads"
 
qc_reads="/mnt/c/Users/Dell/Desktop/yeast/qc_reads"    # for qc results 

## Threads for fastp 
THREADS=4


### Step 1: QC using fastp ###

for R1 in $reads/*_1.fastq.gz; do
     Identify corresponding R2
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    
    # Get sample name
    SAMPLE=$(basename $R1 _1.fastq.gz)
    
    # Output files
    OUT1=$qc_reads/${SAMPLE}_1.clean.fastq.gz
    OUT2=$qc_reads/${SAMPLE}_2.clean.fastq.gz
    HTML=$qc_reads/${SAMPLE}_fastp.html
    JSON=$qc_reads/${SAMPLE}_fastp.json
    
    echo "Running fastp for sample: $SAMPLE"
    
    fastp -i $R1 -I $R2 \
          -o $OUT1 -O $OUT2 \
          -h $HTML -j $JSON \
          -w $THREADS \
          --detect_adapter_for_pe \
          --trim_front1 3 --trim_front2 3 \
          --cut_tail
done

echo "QC complete. Cleaned FASTQ and reports are in $qc_reads"


## Step 2 : Alignment using the  BWA MEM tool ##

# Sample name (assuming single sample)
R1=$qc_reads/*_1.clean.fastq.gz
R2=$qc_reads/*_2.clean.fastq.gz
SAMPLE=$(basename $R1 _1.clean.fastq.gz)

echo "Starting alignment for sample: $SAMPLE"

# 1️ Align with BWA MEM
bwa mem -t 4 $ref $R1 $R2 | \
    samtools view -Sb - > $aligned_reads/${SAMPLE}.bam

# 2️ Sort BAM
samtools sort -@ 1 -m 300M -T /mnt/c/Users/Dell/Desktop/yeast/tmp/${SAMPLE} \       # due to RAM issue 
    $aligned_reads/${SAMPLE}.bam -o $aligned_reads/${SAMPLE}.sorted.bam

# 3️ Index BAM
samtools index $aligned_reads/${SAMPLE}.sorted.bam

echo "Alignment complete. Sorted and indexed BAM is in $aligned_reads"


## Step 3: Variant calling using bcftools  ## 

# Sample BAM (sorted & indexed)

BAM="$aligned_reads/SRR21747980.sorted.bam"

# Output VCF
VCF="$data/SRR21747980.vcf.gz"

# 1️ Generate pileup and call variants
bcftools mpileup -f $ref $BAM -Ou | bcftools call -mv -Oz --ploidy 1 -o $VCF

# 2️ Index the VCF
bcftools index $VCF

echo "Variant calling complete. VCF is at $VCF"

## Step 4 : Variant filtering  ## 

# Input VCF
VCF="/mnt/c/Users/Dell/Desktop/yeast/data/SRR21747980.vcf.gz"

# Output filtered VCF
FILTERED_VCF="/mnt/c/Users/Dell/Desktop/yeast/data/SRR21747980.filtered.vcf.gz"

# Filter by quality and depth
bcftools filter -s LOWQUAL -e 'QUAL<20 || INFO/DP<10' "$VCF" -Oz -o "$FILTERED_VCF"

# Index filtered VCF
bcftools index $FILTERED_VCF


##  SNPS and indels separate ## 

SNP_VCF="$data/SRR21747980.filtered.snps.vcf.gz"
INDEL_VCF="$data/SRR21747980.filtered.indels.vcf.gz"

# SNPs
bcftools view -v snps $data/SRR21747980.filtered.vcf.gz -Oz -o $SNP_VCF
bcftools index $SNP_VCF

# INDELs
bcftools view -v indels $data/SRR21747980.filtered.vcf.gz -Oz -o $INDEL_VCF
bcftools index $INDEL_VCF

echo "SNPs and INDELs separation complete."
