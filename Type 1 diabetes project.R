library(GEOquery)
library(dplyr)
library(readr)

# load the raw count data 
raw_data_1 <- read.table(gzfile("C:/Users/ragha/Desktop/diabetes project/GSE123658_raw_counts_GRCh38.p13_NCBI (2).tsv.gz"), header = FALSE, sep = "\t") 
raw_data_1

# load the annotation data
file_path <- "C:/Users/ragha/Desktop/diabetes project/Human.GRCh38.p13.annot.tsv.gz"
annotation_data_1 <- read_tsv(file_path)

# create a mapping of geneID to symbol (get geneID and sample IFsymbol from annotation data)
geneID_to_symbol_1 <- annotation_data_1 %>% dplyr::select(GeneID, Symbol)
geneID_to_symbol_1

library(tibble)

# renames the columns in raw_data and sets the first row as the column names
colnames(raw_data_1) <- raw_data_1[1, ]
raw_data_1 <- raw_data_1[-1, ]
colnames(raw_data_1)[1] <- "GeneID"
raw_data_1

# Ensure GeneID columns are the same type (character) for the join
raw_data_1$GeneID <- as.character(raw_data_1$GeneID)
geneID_to_symbol_1$GeneID <- as.character(geneID_to_symbol_1$GeneID)

# Merge raw data with annotation(geneID_to_symbol) to get gene symbols
merged_data_1 <- raw_data_1 %>%
  inner_join(geneID_to_symbol_1, by = "GeneID")

# Rearrange columns so that Symbol is the first column
merged_data_1 <- merged_data_1 %>%
  relocate(Symbol, .after = GeneID)

# Remove the GeneID column, if desired so we have the Symbol now
merged_data_1 <- merged_data_1 %>% dplyr::select(-GeneID)
merged_data_1

# Convert count columns to numeric, ignoring any potential warnings
merged_data_1[,-1] <- lapply(merged_data_1[,-1], as.numeric)

# remove duplicates by summing counts for each gene symbol
merged_data_1 <- merged_data_1 %>%
  group_by(Symbol) %>%
  summarize(across(everything(), sum, na.rm = TRUE))
#merged_data now has the raw data with symbol in fist column V

# Removing genes with low counts across all samples
# Determine a minimum count threshold
min_count <- 10
min_samples <- ncol(merged_data_1) / 2

# Filter out genes with low counts
filtered_merged_raw_counts_1 <- merged_data_1[rowSums(merged_data_1 >= min_count) >= min_samples, ]
filtered_merged_raw_counts_1

# rename the symbol column to geneID
colnames(filtered_merged_raw_counts_1)[colnames(filtered_merged_raw_counts_1) == "Symbol"] <- "GeneID"


###


# Load the metadata file 
lines <- readLines(gzfile("C:/Users/ragha/Desktop/diabetes project/GSE123658-GPL18573_series_matrix (4).txt.gz"))
# Extract metadata lines starting with '!Sample_'
metadata_lines <- grep("^!Sample_", lines, value = TRUE)

# Split by tabs and bind as a data frame
metadata_file <- do.call(rbind, strsplit(metadata_lines, "\t"))
metadata_file <- as.data.frame(metadata_file, stringsAsFactors = FALSE)

# View the resulting sample metadata
metadata_file

# Transpose the data frame
metadata_transposed_1 <- t(metadata_file)
metadata_transposed_1 <- as.data.frame(metadata_transposed_1, stringsAsFactors = FALSE)

# Set the first row as column names
colnames(metadata_transposed_1) <- metadata_transposed_1[1, ]
metadata_transposed_1 <- metadata_transposed_1[-1, ]

# Rename columns as needed
colnames(metadata_transposed_1)[10:18] <- c("inventory_patient_id", "sex", "cell_type", "diagnosis", "apoe", "expired_age", "pmi", "braak_score")

# Convert relevant columns to factors variables, this is essential in R when preparing data for certain types of analyses, such as differential expression with DESeq2
metadata_transposed_1$diagnosis <- as.factor(metadata_transposed_1$diagnosis)
if ("batch" %in% colnames(metadata_transposed_1)) {
  metadata_transposed_1$batch <- as.factor(metadata_transposed_1$batch)
}
metadata_transposed_1$sex <- as.factor(metadata_transposed_1$sex)

# Save transformed metadata
write.csv(metadata_transposed_1, "transformed_metadata_file.csv", row.names = FALSE)
metadata_transposed_1

# to see where the file got saved in my computer
getwd()


#now:

# Column names of the count table (excluding GeneID)
colnames(filtered_merged_raw_counts_1)

# Sample identifiers from metadata
colnames(metadata_transposed_1)
metadata_transposed_1[["!Sample_geo_accession"]]

# Remove any unnecessary characters (quotes, spaces, etc.) from metadata sample identifiers.
metadata_transposed_1[["!Sample_geo_accession"]] <- gsub('"', '', metadata_transposed_1[["!Sample_geo_accession"]])
metadata_transposed_1[["!Sample_geo_accession"]] <- trimws(metadata_transposed_1[["!Sample_geo_accession"]])
metadata_transposed_1[["!Sample_geo_accession"]]


# Use double brackets to reference the column
metadata_transposed_1[["!Sample_geo_accession"]] <- gsub("old_pattern", "new_pattern", metadata_transposed_1[["!Sample_geo_accession"]])

# Rename the column in the metadata data frame
colnames(metadata_transposed_1)[colnames(metadata_transposed_1) == "!Sample_geo_accession"] <- "Sample_geo_accession"

# Now you can access it easily using $
metadata_transposed_1$Sample_geo_accession <- gsub("old_pattern", "new_pattern", metadata_transposed_1$Sample_geo_accession)
print(metadata_transposed_1$Sample_geo_accession)

metadata_transposed_1[["!Sample_title"]] <- gsub('"', '', metadata_transposed_1[["!Sample_title"]])
metadata_transposed_1[["!Sample_title"]] <- trimws(metadata_transposed_1[["!Sample_title"]])
metadata_transposed_1[["!Sample_title"]]

colnames(metadata_transposed_1)[colnames(metadata_transposed_1) == "!Sample_title"] <- "Sample_title"

metadata_transposed_2 <- metadata_transposed_1[, c("Sample_title", "Sample_geo_accession")]

####

# Set 'Sample_geo_accession' as row names
rownames(metadata_transposed_2) <- metadata_transposed_2$Sample_geo_accession


# Rename "healthy" to "Healthy" and anything else to "T1D"
if ("Sample_title" %in% colnames(metadata_transposed_2)) {
  metadata_transposed_2$Sample_title <- ifelse(grepl("healthy", metadata_transposed_2$Sample_title, ignore.case = TRUE), 
                              "Healthy", 
                              "T1D")
} else {
  stop("Column 'Sample_title' not found in the data.")
}


# Store GeneID values as row names
filtered_merged_raw_counts_1 <- filtered_merged_raw_counts_1 %>%
  column_to_rownames(var = "GeneID")


# check the dimension to be match between the raw_final and filtered_metadata 
dim(metadata_transposed_2)
dim(filtered_merged_raw_counts_1)


# Extract sample IDs from metadata_transposed_2 (assuming they are in the second column)
metadata_samples <- metadata_transposed_2[[2]]  # Replace '2' if the sample IDs are in a different column

# Extract column names from the filtered merged raw counts (excluding the first column if it's a gene identifier)
row_count_samples <- colnames(filtered_merged_raw_counts_1)[-1]  

# Find common sample IDs
common_samples <- intersect(metadata_samples, row_count_samples)

# Check the name of the first column in the filtered raw counts
first_column_name <- colnames(filtered_merged_raw_counts_1)[1]  # Retrieve the actual name of the first column

# Subset filtered raw counts to include only common samples (and keep the first column for gene IDs)
filtered_counts_matched <- filtered_merged_raw_counts_1[, c(first_column_name, common_samples), drop = FALSE]

# Subset metadata to include only rows with common samples
metadata_matched <- metadata_transposed_2[metadata_transposed_2[[2]] %in% common_samples, ]

# Remove the first column from the filtered_counts_matched data frame
filtered_counts_matched_2 <- filtered_counts_matched[, -1]

# check the dimension to be match between the raw_final and filtered_metadata 
dim(filtered_counts_matched_2)
dim(metadata_matched)


#this is the final preparation 

#filtered_counts_matches_2 and metadata_matched are the files we going to use for analysis 

###################################
#MA plot
#volcano plot
#heatmap
#bar plot
#KEGG
#GO
#pca

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts_matched_2,
  colData = metadata_matched,
  design = ~ Sample_title 
)


# Run DESeq2 normalization and differential expression
dds <- DESeq(dds)

# Results for the condition of interest (e.g., T1D vs Healthy)
res <- results(dds, contrast = c("Sample_title", "T1D", "Healthy"))

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Filter results for significance (adjusted p-value < 0.05)
sig_genes <- res_ordered[which(res_ordered$padj < 0.05), ]

#MA plot 
plotMA(res, ylim = c(-5, 5), main = "MA Plot")



# Filter for significant genes with a log2FoldChange threshold
key_genes <- sig_genes[abs(sig_genes$log2FoldChange) > 1, ]  # Log2 fold change > 1 or < -1
key_genes <- as.data.frame(key_genes)  # Convert to data frame for easier handling

# View the top key genes
head(key_genes)

# Ensure row names in `key_genes` match the gene symbols
key_genes$GeneID <- rownames(key_genes)

# Merge with your annotation file to include gene symbols
key_genes_annotated <- merge(key_genes, geneID_to_symbol_1, by.x = "GeneID", by.y = "GeneID")
key_genes_annotated <- key_genes_annotated[order(key_genes_annotated$padj), ]  # Sort by adjusted p-value

# View the annotated key genes
head(key_genes_annotated)



# Enhanced MA plot with highlighted significant genes
plotMA(res, ylim = c(-5, 5), main = "Enhanced MA Plot")
points(
  x = log10(key_genes$baseMean), 
  y = key_genes$log2FoldChange, 
  col = "red", 
  pch = 20
)


# Count the number of upregulated and downregulated genes
upregulated <- sum(key_genes$log2FoldChange > 1)
downregulated <- sum(key_genes$log2FoldChange < -1)

cat("Number of upregulated genes:", upregulated, "\n")
cat("Number of downregulated genes:", downregulated, "\n")

#number of upregulated genes are 46
#number of downregulated genes are 30  

# Create a data frame for the counts
regulation_summary <- data.frame(
  Regulation = c("Upregulated", "Downregulated"),
  Count = c(46, 30)
)

# Load ggplot2 for visualization
library(ggplot2)

# Create the bar chart to visualize the upregulated genes and downregulated ones
ggplot(regulation_summary, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "Number of Differentially Expressed Genes",
       x = "Regulation Type",
       y = "Count") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  )



# volcano plot 
library(ggplot2)

# Create a data frame for plotting
res_volcano <- as.data.frame(res)
res_volcano$GeneID <- rownames(res_volcano)  # Add GeneID for labeling

# Add significance and fold change thresholds
res_volcano$Significance <- ifelse(
  res_volcano$padj < 0.05 & abs(res_volcano$log2FoldChange) > 1, 
  ifelse(res_volcano$log2FoldChange > 0, "Upregulated", "Downregulated"), 
  "Not Significant"
)

# Plot
ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  labs(
    title = "Volcano Plot", 
    x = "Log2 Fold Change", 
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

############

#heatmap 

# Extract normalized counts for significant genes
normalized_counts <- counts(dds, normalized = TRUE)
heatmap_data <- normalized_counts[rownames(normalized_counts) %in% rownames(sig_genes), ]

# Scale the data for visualization
heatmap_data_scaled <- t(scale(t(heatmap_data)))  # Z-score normalization

# Create heatmap
pheatmap(
  heatmap_data_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = metadata_matched,  # Add sample annotations if available
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  fontsize = 8,
  main = "Heatmap of Significant Genes"
)



###################

#PCA
# Perform PCA
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "Sample_title", returnData = TRUE)

# Custom PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample_title)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")

#######

#KEGG and GO Enrichment Analysis

library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs (required for enrichment analysis)
sig_genes_ids <- bitr(
  rownames(sig_genes),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# KEGG pathway enrichment
kegg_results <- enrichKEGG(
  gene = sig_genes_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# GO enrichment (biological processes)
go_results <- enrichGO(
  gene = sig_genes_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological processes
  pvalueCutoff = 0.05
)

# Visualize KEGG results
dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")

# Visualize GO results
dotplot(go_results, showCategory = 10, title = "GO Biological Processes")


##############

#Principal Component Analysis (PCA)
#To visualize the overall variation in the samples:

# Transform the data using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
plotPCA(vsd, intgroup = "Sample_title") +
  ggtitle("PCA of Samples")

##########


#####

# Ensure ggplot2 is loaded
library(ggplot2)

# Extract top 10 significant genes
top10_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
top10_genes_df <- as.data.frame(top10_genes)
top10_genes_df$GeneID <- rownames(top10_genes_df)

# Create a bar chart of log2 fold changes for the top 10 genes
ggplot(top10_genes_df, aes(x = reorder(GeneID, log2FoldChange), y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for better readability
  theme_minimal() +
  labs(
    title = "Top 10 Significant Genes",
    x = "Gene ID",
    y = "Log2 Fold Change"
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Extract top 10 significant genes by adjusted p-value
top10_genes <- head(sig_genes[order(sig_genes$padj), ], 10)

# Add gene symbols as a column (if they are row names)
top10_genes_df <- as.data.frame(top10_genes)
top10_genes_df$GeneID <- rownames(top10_genes_df)

# View the top 10 genes
print(top10_genes_df)

# Save the results to a CSV file
write.csv(top10_genes_df, "Top10_Significant_Genes.csv", row.names = FALSE)


##########








