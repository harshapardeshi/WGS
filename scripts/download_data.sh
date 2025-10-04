#!/bin/bash 

## Before running, make sure the following directories exist:
## for reads:-  /mnt/c/Users/Dell/Desktop/yeast/reads
## for ref  :-  /mnt/c/Users/Dell/Desktop/yeast/ref
## You can create them using:
## mkdir -p /mnt/c/Users/Dell/Desktop/yeast/reads
## mkdir -p /mnt/c/Users/Dell/Desktop/yeast/ref

## download the data ##

wget -P /mnt/c/Users/Dell/Desktop/yeast/reads/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR217/080/SRR21747980/SRR21747980_1.fastq.gz
wget -P /mnt/c/Users/Dell/Desktop/yeast/reads/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR217/080/SRR21747980/SRR21747980_2.fastq.gz

echo "Run the Prep files..."

# Download reference
wget -P /mnt/c/Users/Dell/Desktop/yeast/ref/ \
ftp://ftp.ensemblgenomes.org/pub/fungi/release-62/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

# Unzip
gunzip -f /mnt/c/Users/Dell/Desktop/yeast/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

# index the ref using bwa   
bwa index /mnt/c/Users/Dell/Desktop/yeast/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

# samtools FASTA index
samtools faidx  /mnt/c/Users/Dell/Desktop/yeast/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

