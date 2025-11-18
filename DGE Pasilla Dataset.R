# Install packages (run once)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "pasilla", "pheatmap",
                       "EnhancedVolcano"))


# ==============================
# 2. Load libraries
# ==============================
library(DESeq2)
library(pasilla)
library(pheatmap)
library(EnhancedVolcano)


# ==============================
# 3. Load pasilla count data
# ==============================
datafile <- system.file("extdata/pasilla_gene_counts.tsv",
                        package="pasilla")
raw_counts <- read.table(datafile, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE)


# Keep only sample columns (7 samples)
countData <- raw_counts[, c("untreated1", "untreated2",
                            "untreated3", "untreated4",
                            "treated1", "treated2", "treated3")]
# Convert to numeric matrix safely
# Should be numeric


# ==============================
# 4. Load and align metadata
# ==============================
annotation <-
  system.file("extdata/pasilla_sample_annotation.csv",
              package="pasilla")
colData <- read.csv(annotation, row.names=1)

# Rename rownames to match countData columns
rownames(colData) <-
  c("untreated1","untreated2","untreated3","untreated4",
    "treated1","treated2","treated3")


# Reorder rows to match countData columns
colData <- colData[colnames(countData), ]


# Verify alignment
all(colnames(countData) == rownames(colData)) 
# Should return TRUE

colData <- colData["condition"]
colData$condition <- factor(colData$condition)


# ==============================
# 5. Create DESeq2 dataset
# ==============================
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]


# ==============================
# 6. Run DESeq2
# ==============================
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]
head(resOrdered)


# Save results
write.csv(as.data.frame(resOrdered),
          "DESeq2_pasilla_results.csv")


# ==============================
# 7. Visualization
# ==============================
# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)



# PCA plot
plotPCA(vsd, intgroup="condition")
# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample-to-Sample Distances")


# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Pasilla Knockdown vs Control')
