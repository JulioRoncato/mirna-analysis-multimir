# Load required libraries
library(ggplot2)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(VennDiagram)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(png)
library(grid)

# Function to create figures and save them as PNG
create_mirna_figures <- function(){
  # Assuming data is loaded or passed as parameters

  # 1. Generate Volcano Plot
  volcano_data <- data.frame(logFC = rnorm(100), pval = runif(100))  # Replace with actual miRNA differential expression data
  p1 <- ggplot(volcano_data, aes(x = logFC, y = -log10(pval))) +
    geom_point() +
    theme_minimal() +
    ggtitle("Volcano Plot of miRNA Differential Expression")
  ggsave("figures/volcano_plot.png", plot = p1, height = 5, width = 7, dpi = 300)

  # 2. Create miRNA-target network visualization
  # Example data for validated and predicted targets
  edges <- data.frame(from = c("miRNA1", "miRNA2"), to = c("target1", "target2"))  # Replace with actual data
  g <- graph_from_data_frame(edges)
  p2 <- ggraph(g) +
    geom_edge_link() +
    geom_node_point() +
    theme_minimal() +
    ggtitle("miRNA-Target Network Visualization")
  ggsave("figures/miRNA_target_network.png", plot = p2, height = 5, width = 7, dpi = 300)

  # 3. Generate GO enrichment dot plot
  go_results <- data.frame(term = c("Term1", "Term2"), pvalue = c(0.01, 0.02))  # Replace with actual GO results
  p3 <- ggplot(go_results, aes(x = term, y = -log10(pvalue))) +
    geom_point() +
    theme_minimal() +
    ggtitle("GO Enrichment Dot Plot")
  ggsave("figures/go_enrichment_dot_plot.png", plot = p3, height = 5, width = 7, dpi = 300)

  # 4. Create Venn Diagram
  valid_targets <- c("target1", "target2")
  pred_targets <- c("target2", "target3")  # Replace with actual target lists
  venn.plot <- venn.diagram(
    x = list(Validated = valid_targets, Predicted = pred_targets), 
    filename = NULL,
    output = TRUE
  )
  png("figures/venn_diagram.png", width = 7, height = 5, units = "in", res = 300)
  grid.draw(venn.plot)
  dev.off()

  # 5. Generate Heatmap
  heatmap_data <- matrix(rnorm(100), nrow = 10)  # Replace with actual heatmap data
  p5 <- Heatmap(heatmap_data)
  png("figures/heatmap.png", width = 7, height = 5, units = "in", res = 300)
  draw(p5)
  dev.off()
}

# Call the function to create the figures
create_mirna_figures()