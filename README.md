# 🔬 microRNA (miRNA) Analysis - Disease Biomarkers

## 🧬 Biological Background

microRNAs (miRNAs) are small non-coding RNAs (~22 nt) that regulate gene expression post-transcriptionally by binding to mRNA 3'UTRs. Dysregulated miRNA expression is associated with:
- **Cancer progression** (oncomiRs, tumor suppressors)
- **Cardiovascular disease**
- **Neurodegenerative disease**
- **Inflammatory conditions**

This project characterizes miRNA expression profiles and identifies disease-associated miRNAs and their regulatory networks using **validated and predicted targets**.

---

## 🧪 Experimental Design

| Analysis Type | Goal |
|---|---|
| Small RNA-seq | Quantify miRNA expression across samples |
| Quality Control | Assess library complexity and artifact presence |
| Normalization | Remove technical variability |
| Differential Expression | Identify dysregulated miRNAs |
| Target Retrieval | Obtain validated and predicted miRNA targets |
| Target Validation | Filter high-confidence interactions |
| Gene Enrichment | Pathway analysis of miRNA-regulated genes |

---

## ⚙️ Bioinformatics Workflow

### 1️⃣ Read Alignment & Quantification
- Adapter trimming with **Trimmomatic** or **Cutadapt**
- Alignment to miRNA reference (miRBase) using **bowtie2** or **HISAT2**
- Quantification of read counts per miRNA with **featureCounts** or **miRDeep2**

### 2️⃣ Quality Control
- Read length distribution
- Adapter contamination
- rRNA/tRNA contamination
- Library complexity metrics

### 3️⃣ Normalization
- Library size normalization (TMM, DESeq2)
- Batch effect removal (ComBat, SVA)
- Variance stabilization

### 4️⃣ Differential Expression Analysis
- **DESeq2** or **edgeR** for differential miRNA expression
- FDR-corrected p-values (adj. p < 0.05)
- Log2 fold-change filtering (|log2FC| > 1)

### 5️⃣ miRNA Target Retrieval & Validation

**This is where `multiMiR` becomes essential:**

#### 5A. Target Database Integration with **multiMiR**

The `multiMiR` package integrates multiple target prediction and validation databases:

```r
# Install and load multiMiR
# BiocManager::install("multiMiR")
library(multiMiR)

# Get validated targets from multiple databases
# Validated interactions (experimentally confirmed)
validated_targets <- get_multimir(
  org = "hsa",
  mirna = c("hsa-miR-155", "hsa-miR-146a"),
  table = "validated"
)

# Predicted targets from multiple algorithms
predicted_targets <- get_multimir(
  org = "hsa",
  mirna = c("hsa-miR-155", "hsa-miR-146a"),
  table = "predicted"
)

# Combined: validated + predicted
all_targets <- get_multimir(
  org = "hsa",
  mirna = c("hsa-miR-155", "hsa-miR-146a"),
  table = "all"
)
```

#### 5B. Integrated Target Datasets

`multiMiR` retrieves targets from:

| Target Type | Databases Included |
|---|---|
| **Validated Targets** | miRTarBase, miRecords, MIRT |
| **Predicted Targets** | TargetScan, DIANA-microT, PicTar, PITA |
| **Metadata** | Target genes, binding sites, publications |

#### 5C. Target Filtering Strategy

```r
# Filter validated targets (high confidence)
validated_df <- validated_targets@data %>%
  filter(experiment_support > 0)  # At least 1 experimental confirmation

# For predicted targets, apply stringency criteria
predicted_df <- predicted_targets@data %>%
  filter(predicted_score > 0.5)  # Example: algorithm confidence score

# Combine: prioritize validated, supplement with high-confidence predicted
high_confidence_targets <- rbind(
  validated_df %>% mutate(confidence = "validated"),
  predicted_df %>% mutate(confidence = "predicted")
)
```

### 6️⃣ Target Gene Enrichment Analysis

From the gene set of validated + predicted targets:

#### 6A. Gene Ontology (GO) Enrichment
```r
library(clusterProfiler)

# Extract unique genes from targets
target_genes <- unique(high_confidence_targets$target_gene)

# GO enrichment
go_enrichment <- enrichGO(
  gene = target_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
```

#### 6B. KEGG Pathway Enrichment
```r
# Convert gene symbols to Entrez IDs
gene_ids <- bitr(target_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# KEGG enrichment
kegg_enrichment <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
```

#### 6C. Disease-Specific Enrichment
```r
# Using DisGeNET or similar disease databases
# Map target genes to disease associations
disease_enrichment <- enrichDGN(
  gene = gene_ids$ENTREZID,
  pvalueCutoff = 0.05
)
```

### 7️⃣ Regulatory Network Construction

```r
library(igraph)
library(ggraph)

# Build miRNA-gene network
# Nodes: dysregulated miRNAs + target genes
# Edges: miRNA-target interactions (validated vs predicted)

network <- graph_from_data_frame(
  d = high_confidence_targets %>% 
      select(miRNA = mirna, Gene = target_gene, confidence),
  directed = TRUE
)

# Visualize with edge coloring by confidence
plot(network, 
     edge.color = ifelse(E(network)$confidence == "validated", 
                         "red", "gray"),
     edge.width = 2,
     vertex.size = 10)
```

---

## 📊 Data Visualization

### miRNA Expression Heatmap
![miRNA Heatmap](figures/mirna_heatmap.png)

### Differential Expression Volcano Plot
![Volcano Plot](figures/mirna_volcano.png)

### multiMiR Target Integration
**Venn Diagram**: Validated vs Predicted targets overlap
![Target Venn](figures/mirna_target_venn.png)

### Target Gene Enrichment
![GO Enrichment](figures/go_enrichment_dotplot.png)

### Regulatory Network
**Network Visualization**: miRNA → Gene connections (colored by target validation status)
![Network](figures/mirna_target_network.png)

---

## 🔐 Data Availability

Small RNA-seq data from public repositories:
- **GEO** (Gene Expression Omnibus)
- **SRA** (Sequence Read Archive)
- **TCGA** (The Cancer Genome Atlas)

Target databases accessed via `multiMiR`:
- miRTarBase, miRecords, MIRT (validated)
- TargetScan, DIANA-microT, PicTar, PITA (predicted)

---

## 🛠️ Tools and Software Used

### Alignment & Quantification
- **Cutadapt** - Adapter trimming
- **bowtie2** - Short-read alignment
- **miRDeep2** - miRNA discovery and quantification
- **featureCounts** - Read counting

### R / Bioconductor
- **DESeq2** - Differential expression analysis
- **edgeR** - Alternative differential expression
- **multiMiR** - ⭐ **Integrated miRNA target database retrieval and management**
  - Retrieves validated targets from miRTarBase, miRecords, MIRT
  - Retrieves predicted targets from TargetScan, DIANA-microT, PicTar, PITA
  - Filters by confidence, experiment type, and prediction algorithms
  - Returns curated datasets with metadata
- **Biostrings** - Sequence manipulation
- **clusterProfiler** - Functional enrichment (GO, KEGG)
- **org.Hs.eg.db** - Human gene annotations
- **igraph** / **ggraph** - Network analysis and visualization
- **ggplot2** - Data visualization

### Databases (via multiMiR)
- **miRBase** (http://mirbase.org) - miRNA sequences and annotations
- **miRTarBase** (http://mirtarbase.cuhk.edu.cn) - Validated targets
- **TargetScan** (http://www.targetscan.org) - Predicted targets
- **DIANA-microT** - microT scores
- **PicTar**, **PITA** - Additional prediction algorithms

---

## 📌 Clinical Relevance

This project demonstrates:

- **Comprehensive target identification**: Integration of validated + predicted targets provides complete regulatory landscape
- **High-confidence biomarker discovery**: Prioritization of experimentally confirmed interactions
- **Mechanistic insights**: Understanding miRNA-mediated gene regulation pathways
- **Therapeutic targets**: 
  - miRNAs as drug targets (anti-miRs, miRNA inhibitors)
  - Target genes as downstream therapeutic opportunities
- **Disease stratification**: miRNA signatures correlate with disease subtypes and outcomes
- **Precision medicine**: Patient stratification based on miRNA-target-pathway networks

### Example Clinical Application
```
Dysregulated miR-155-5p
    ↓
multiMiR retrieval → 127 validated targets + 450 predicted targets
    ↓
Filter to high-confidence → 89 targets with multiple confirmations
    ↓
GO/KEGG enrichment → Pathways: NF-κB signaling, inflammatory response, immune cell activation
    ↓
Clinical interpretation → miR-155 drives inflammation in disease X
    ↓
Therapeutic opportunity → Targeting miR-155 or its validated targets (e.g., SOCS1, SHIP1)
```

---

## ⚠️ Considerations & Best Practices

### multiMiR-Specific
- **Database Versions**: Target databases update regularly; specify `multiMiR` version and database snapshot date
- **Cross-Species**: Use appropriate organism code (`org = "hsa"` for human, `"mmu"` for mouse)
- **Target Overlap**: Same target may appear in multiple databases with different confidence levels
- **Redundancy**: Paralogous miRNAs may share targets; filter carefully

### General
- **Reference Database**: miRBase versions differ; specify version used
- **Cross-Reactive Reads**: Some reads may map to multiple miRNAs
- **Tissue/Cell-Type Specificity**: miRNA profiles vary dramatically by tissue
- **Validation**: Computational predictions (even validated ones) require experimental validation (RT-qPCR, Luciferase assays)
- **Circulating vs Tissue miRNAs**: Different extraction and normalization strategies
- **Multiple Testing Correction**: Apply FDR adjustment for enrichment analysis

---

## 📚 Key References for multiMiR

- **multiMiR Package**: https://bioconductor.org/packages/multiMiR
- Ru et al. (2014): "The multiMiR R package and database: integration of microRNA–target interactions along with their disease and drug associations"
- miRTarBase: https://mirtarbase.cuhk.edu.cn
- TargetScan: http://www.targetscan.org
