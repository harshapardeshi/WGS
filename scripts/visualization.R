################################################################################
# Title: Yeast SNP Analysis and Functional Enrichment
# Author: Harsha
# Date: 2025-10-04
# Description: Clean VEP TXT output, visualize SNP distribution, identify
#              top genes, and perform KEGG enrichment in Saccharomyces cerevisiae
################################################################################

# =========================
# 0. Load data
# =========================
# If your TXT is tab-delimited
vep_data <- read.table("annotated.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# If your TXT is space-delimited, use:
# vep_data <- read.table("annotated.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

# Check the first few rows
head(vep_data)

# Check column names
colnames(vep_data)

# Save as TSV for record
write.table(vep_data, "annotated_vep.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Read it back to check
vep_tsv <- read.table("annotated_vep.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(vep_tsv)
colnames(vep_tsv)

# =========================
# 1. Load libraries
# =========================
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("org.Sc.sgd.db")
library(org.Sc.sgd.db)  # Yeast gene annotations

# =========================
# 2. Clean and format data
# =========================
vep_clean <- data.frame(
  Location    = vep_tsv$`I.63.63`,
  Allele      = vep_tsv$C,
  Consequence = vep_tsv$`upstream_gene_variant`,
  Gene        = vep_tsv$`YAL067W.A`,
  stringsAsFactors = FALSE
)

head(vep_clean)

# Split Location into Chromosome and Position
vep_clean <- vep_clean %>%
  separate(Location, into = c("Chr", "Pos"), sep = "[:-]", convert = TRUE)

head(vep_clean)

# =========================
# 3. SNP Summary
# =========================
snp_count <- vep_clean %>% 
  group_by(Chr) %>% 
  summarise(count = n())
print(snp_count)

consequence_count <- vep_clean %>% 
  group_by(Consequence) %>% 
  summarise(count = n())
print(consequence_count)

# =========================
# 4. Visualization
# =========================
# SNP count per chromosome
ggplot(snp_count, aes(x = Chr, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "SNP Count per Chromosome in Yeast",
       x = "Chromosome",
       y = "Number of SNPs")

# SNP consequences
ggplot(consequence_count, aes(x = reorder(Consequence, -count), y = count)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  coord_flip() +
  labs(title = "SNP Consequences in Yeast",
       x = "Consequence Type",
       y = "Count")

# =========================
# 5. Top 20 genes by SNP count
# =========================
top_genes <- vep_clean %>%
  group_by(Gene) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  head(20)

# Plot top 20 genes
ggplot(top_genes, aes(x = reorder(Gene, count), y = count)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  labs(title = "Top 20 Genes by SNP Count",
       x = "Gene",
       y = "Number of SNPs")

# =========================
# 6. Functional Enrichment (KEGG)
# =========================
gene_list <- top_genes$Gene
gene_list

kegg_enrich <- enrichKEGG(
  gene         = gene_list,
  organism     = "sce",   # Saccharomyces cerevisiae
  keyType      = "kegg",
  pvalueCutoff = 0.05
)

# View results
head(kegg_enrich)

# KEGG dotplot
dotplot(kegg_enrich, showCategory = 10, title = "KEGG Pathways")

# Save KEGG enrichment results
write.csv(as.data.frame(kegg_enrich), "KEGG_enrichment_top20_genes.csv", row.names = FALSE)

################################################################################
# End of Script
################################################################################
