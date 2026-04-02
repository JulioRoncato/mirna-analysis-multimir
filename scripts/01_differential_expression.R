# Load necessary libraries
library(DESeq2)
library(edgeR)

# Read input data
# Assuming raw counts are stored in a CSV file
countData <- read.csv("path/to/counts.csv", row.names = 1)
colData <- read.csv("path/to/colData.csv")

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition) # 'condition' should be the column in colData

dds <- DESeq(dds)
res_deseq <- results(dds)

# Save DESeq2 results
write.csv(as.data.frame(res_deseq), file = "DESeq2_results.csv")

# edgeR analysis
group <- factor(colData$condition) # 'condition' should be the column in colData
y <- DGEList(counts = countData, group = group)
y <- calcNormFactors(y)
design <- model.matrix(~ group)

# Estimate dispersion
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

# Save edgeR results
topTags <- topTags(lrt, n = Inf)
write.csv(as.data.frame(topTags), file = "edgeR_results.csv")
