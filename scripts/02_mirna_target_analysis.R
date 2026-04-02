# Load required libraries\nlibrary(multiMiR)
library(clusterProfiler)\nlibrary(org.Hs.eg.db)  # Use the appropriate organism package
# Function to retrieve targets
retrieve_targets <- function(mirna_ids) {targets <- multiMiR::get_multiMir(miRNA = mirna_ids, species = 'hsa',  # Change as needed for other species  
                                                                           validation = TRUE)
                                                                           return(targets)}
# Filter targets by confidence level
filter_targets <- function(targets_df, confidence_threshold) {filtered_targets <- targets_df[targets_df$confidence >= confidence_threshold, ]  
                                                              return(filtered_targets)}

# Perform gene enrichment analysis
perform_enrichment <- function(genes) {  
                                       go_results <- enrichGO(
                                         gene = genes,\n    OrgDb = org.Hs.eg.db,\n    keyType = 'ENTREZID',
                                         pAdjustMethod = 'BH',\n    qvalueCutoff = 0.05\n  )
  
  kegg_results <- enrichKEGG(
    gene = genes,
    organism = 'hsa',
    pvalueCutoff = 0.05
  )
  
  list(go = go_results, kegg = kegg_results)
}

# Example Usage
mirna_ids <- c("miR-21", "miR-155")  

# Example miRNA IDs
confidence_threshold <- 0.5

# Retrieve and filter targets
targets <- retrieve_targets(mirna_ids)
filtered_targets <- filter_targets(targets$data, confidence_threshold)


# Gene enrichment analysis
enrichment_results <- perform_enrichment(filtered_targets$entrezGene)

# Output results
print(enrichment_results)
