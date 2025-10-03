# 
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
# Load necessary libraries
library(data.table)
library(dplyr)
library(limma)
library(pheatmap)
library(ggplot2)
library(AnnotationDbi)  # Annotation package
library(org.Hs.eg.db)   # Human gene annotation package

# Load gene expression data
expr_data <- fread("/N/project/SingleCell_Image/Ronak/anmerge_task/gene_exp.csv")

# Load diagnosis data
pheno_data <- fread("/N/project/SingleCell_Image/Ronak/anmerge_task/combined_df.csv")

# Keep only CN and AD samples
pheno_data <- pheno_data %>% filter(Final_Diagnosis %in% c("CTL", "AD"))

# Extract sample IDs from phenotype data
sample_ids <- pheno_data$gexp_id

# Find common samples
common_samples <- intersect(colnames(expr_data)[-1], sample_ids)

# Filter expression data to keep only common samples
expr_data_filtered <- expr_data[, c("V1", common_samples), with=FALSE]

# Ensure pheno_data also only contains common samples
pheno_data_filtered <- pheno_data %>% filter(gexp_id %in% common_samples)

# Reorder the expression data columns to match the order in phenotype data
expr_data_ordered <- expr_data_filtered[, c("V1", pheno_data_filtered$gexp_id), with=FALSE]

# Map probe IDs to gene symbols using org.Hs.eg.db
probe_ids <- expr_data_ordered$V1
gene_symbols <- mapIds(org.Hs.eg.db, keys = probe_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace probe IDs with gene symbols
expr_data_ordered$V1 <- gene_symbols

# Remove rows with NA gene symbols
expr_data_ordered <- expr_data_ordered[!is.na(expr_data_ordered$V1), ]

# Ensure gene expression data is in matrix form for limma
expr_matrix <- as.matrix(expr_data_ordered[, -1, with=FALSE])
rownames(expr_matrix) <- expr_data_ordered$V1

# Create the design matrix
design <- model.matrix(~0 + factor(pheno_data_filtered$Final_Diagnosis))
colnames(design) <- levels(factor(pheno_data_filtered$Final_Diagnosis))

### Batch effect check
pca <- prcomp(t(expr_matrix), scale.=TRUE)
pca_data <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], Diagnosis = pheno_data_filtered$Final_Diagnosis)

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Diagnosis)) +
  geom_point() +
  labs(title = "PCA Plot", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Perform voom transformation (variance stabilizing transformation)
v <- voom(expr_matrix, design)

# Fit the linear model
fit <- lmFit(v, design)

# Create contrast matrix to compare AD vs CTL
contrast_matrix <- makeContrasts(
  ADvsCTL = AD - CTL,
  levels=design
)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get the results for the contrasts
ADvsCTL <- topTable(fit2, coef="ADvsCTL", number=Inf, p.value=1, sort.by="p")

# Visualize the results

# Heatmap of top differentially expressed genes
top_genes <- ADvsCTL %>% top_n(50, wt=-adj.P.Val)
top_gene_symbols <- top_genes$ID

# Subset the expression matrix to these top genes
heatmap_data <- expr_matrix[rownames(expr_matrix) %in% top_gene_symbols, ]

# Create the heatmap if enough data is present
if (nrow(heatmap_data) >= 2 && ncol(heatmap_data) >= 2) {
  pheatmap(heatmap_data, annotation_col=data.frame(Diagnosis=pheno_data_filtered$Final_Diagnosis),
           show_rownames=TRUE, show_colnames=FALSE)
} else {
  cat("Not enough data to generate a heatmap.\n")
}

# Create a volcano plot
ggplot(ADvsCTL, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=adj.P.Val < 0.05)) +
  theme_minimal() +
  labs(title="Volcano Plot", x="Log Fold Change", y="-log10 Adjusted P-Value") +
  scale_color_manual(values=c("red", "black"))

