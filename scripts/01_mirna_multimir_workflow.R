# miRNA Analysis Script Using multiMiR Package

# Load necessary libraries
library(DESeq2)
library(multiMiR)
library(ggplot2)
library(clusterProfiler)
library(igraph)
library(org.Hs.eg.db)

# Load your count data and metadata
count_data <- read.csv("path_to_your_count_data.csv", row.names = 1)
coldata <- read.csv("path_to_your_metadata.csv", row.names = 1)

# Differential Expression Analysis
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
diff_res <- results(dds)

# Filter results for significant miRNAs
sig_miRNAs <- subset(diff_res, padj < 0.05)

# Target Retrieval with multiMiR
validated_targets <- getMultiMiR(sig_miRNAs$miRNA, species = "Homo sapiens", validatedOnly = TRUE)
predicted_targets <- getMultiMiR(sig_miRNAs$miRNA, species = "Homo sapiens", validatedOnly = FALSE)

# Combine targets
all_targets <- rbind(validated_targets, predicted_targets)

# Filter targets based on user-defined criteria (example: score threshold)
targets_filtered <- all_targets[all_targets$miRNA_score > 0.5, ]

# GO/KEGG Enrichment Analysis
# Convert miRNA targets to gene symbols for enrichment analysis
gene_symbols <- unique(targets_filtered$target)
enrich_go <- enrichGO(gene = gene_symbols, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP', pAdjustMethod = 'BH')
enrich_kegg <- enrichKEGG(gene = gene_symbols, organism = "hsa", pvalueCutoff = 0.05)

# Network Construction
network_data <- data.frame(targets_filtered$miRNA, targets_filtered$target)
network_graph <- graph_from_data_frame(network_data)

# Plot Network
plot(network_graph)

# Generate data for visualization
# Save results
write.csv(as.data.frame(sig_miRNAs), file = "differential_expression_results.csv")
write.csv(as.data.frame(targets_filtered), file = "filtered_targets.csv")
write.csv(as.data.frame(enrich_go), file = "go_enrichment_results.csv")
write.csv(as.data.frame(enrich_kegg), file = "kegg_enrichment_results.csv")

# The End
