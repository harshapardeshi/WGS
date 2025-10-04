# Yeast Whole Genome Sequencing (WGS) Project

Omics Sequencing of *Saccharomyces cerevisiae* Strain with Improved
Capacity for Ethanol Production.

Link:- [MDPI Article](https://www.mdpi.com/2311-5637/9/5/483)

summary This project analyzes whole genome sequencing (WGS) data of
*Saccharomyces cerevisiae* for variant discovery and functional
analysis. The key findings include:

-   **Variant Discovery:** Identified SNPs and INDELs across the genome.
    Filtered variants with `QUAL ≥ 20` and `DP ≥ 10`.\
-   **SNP Consequences:** Most SNPs were located in coding regions, with
    a significant proportion causing missense changes (see
    `snp_consequences.pdf`).\
-   **Chromosomal Distribution:** Variants were unevenly distributed
    across chromosomes, with chromosome I and XII showing higher SNP
    counts (see `snp_count_per_chromosome.pdf`).\
-   **Top Genes by Variants:** Top 20 genes with the highest number of
    SNPs were identified, indicating potential targets for functional
    studies (`top20_genes_by_snp_count.pdf`).\
-   **Pathway Analysis:** KEGG enrichment revealed key metabolic and
    stress-response pathways enriched for genes harboring variants
    (`kegg_pathway.pdf`).

Data acquisition:-
[PRJNA885247](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA885247)

Workflow Data acquition The raw sequencing reads (paired-end) were
downloaded from the European Nucleotide Archive (ENA) using BioProject
ID PRJNA885247. The reference genome and annotation files for
Saccharomyces cerevisiae were obtained from Ensembl Fungi.

QC Quality control (QC) of the FASTQ files was performed using fastp.
The results indicated that the reads were of good quality overall,
making them suitable for downstream analyses.

Alignment The reads were aligned to the reference genome using BWA-MEM.
The resulting alignments were converted to BAM format and subsequently
sorted and indexed using SAMtools.

Variant calling / filtering The aligned reads were used for variant
calling to identify SNPs and INDELs using bcftools. The resulting
variants were stored in a VCF file, containing detailed information
about each variant, including its position, quality, and depth. The VCF
file was then indexed with bcftools for efficient access. High-quality
variants were retained by applying filtering thresholds of QUAL ≥ 20 and
DP ≥ 10. After filtering, SNPs and INDELs were separated into distinct
files for downstream analyses.

Variant annotation Variant annotation was performed using the Ensembl
VEP (Variant Effect Predictor) online tool. The filtered SNP VCF file
was uploaded to VEP, which provided functional consequences for each
variant, including predicted impacts on genes, transcripts, and
proteins. The annotated results were downloaded and saved as a text file
for downstream analysis and visualization.

Visualization The filtered and annotated variants were visualized to
summarize key findings. QC reports from fastp provided overall read
quality, while variant data were used to generate plots such as SNP
consequences, SNP distribution per chromosome, top 20 genes by SNP
count, and KEGG pathway enrichment. These visualizations highlight the
genomic variation and functional impact of variants in Saccharomyces
cerevisiae.
