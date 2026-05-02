
# Bloodborne Biomarkers: Differential RNA Expression in T1D Patients vs. Controls

A bioinformatics project identifying differentially expressed genes (DEGs) in whole blood RNA-seq data from Type 1 Diabetes (T1D) patients versus healthy controls, using the DESeq2 pipeline in R.

> **Author:** Raghad Malkawi  
> **Course project** - Indiana University Luddy School of Informatics, Computation and Engineering, taught by Dr. Juexin Wang, PhD

---

## Background

Type 1 Diabetes Mellitus (T1DM) is a chronic autoimmune disease in which T-cells destroy the pancreatic beta cells responsible for insulin production, leading to insulin deficiency and hyperglycemia. With over 537 million people affected by Diabetes Mellitus worldwide as of 2021, and T1D known to have a significant genetic component, investigating transcriptomic differences between T1D patients and healthy individuals is a meaningful step toward identifying novel biomarkers and better understanding the disease's molecular basis.

---

## Dataset

- **Source:** [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658) — Accession: **GSE123658**
- **Sample type:** Whole blood RNA-seq
- **Samples used:** 46 out of 82 available samples — 32 T1D patients and 14 healthy controls
- **Files needed:**
  - `GSE123658_raw_counts_GRCh38.p13_NCBI.tsv.gz` — Raw gene counts
  - `Human.GRCh38.p13.annot.tsv.gz` — Gene annotation (GeneID → Symbol mapping)
  - `GSE123658-GPL18573_series_matrix.txt.gz` — Sample metadata

All three files can be downloaded directly from the GEO accession page linked above.

---

## Requirements

**R packages:**

| Package | Source | Purpose |
|---|---|---|
| `DESeq2` | Bioconductor | Differential expression analysis |
| `GEOquery` | Bioconductor | GEO data access |
| `clusterProfiler` | Bioconductor | KEGG and GO pathway enrichment |
| `org.Hs.eg.db` | Bioconductor | Human gene ID mapping |
| `pheatmap` | CRAN | Heatmap visualization |
| `RColorBrewer` | CRAN | Color palettes for plots |
| `ggplot2` | CRAN | General plotting |
| `dplyr` | CRAN | Data wrangling |
| `tibble` | CRAN | Data frame utilities |
| `readr` | CRAN | File reading |

Install Bioconductor packages with:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "GEOquery", "clusterProfiler", "org.Hs.eg.db"))
```

Install CRAN packages with:
```r
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tibble", "readr"))
```

---

## Project Structure

```
Diabetes-type-1-project/
├── Type 1 diabetes project.R             # Main analysis script
├── README.md                             # This file
└── outputs/                              # Generated output files (after running)
    ├── transformed_metadata_file.csv
    ├── Top10_Significant_Genes.csv
    └── plots/
        ├── MA_plot.png
        ├── enhanced_MA_plot.png
        ├── volcano_plot.png
        ├── heatmap.png
        ├── bar_plot_DEGs.png
        ├── top10_genes_bar.png
        ├── KEGG_dotplot.png
        ├── GO_dotplot.png
        └── PCA_plot.png
```

---

## Pipeline Overview

### 1. Data Loading & Preprocessing
- Load raw count data and gene annotation file
- Map GeneID to gene symbols via inner join
- Sum counts for duplicate gene symbols
- Filter out low-expression genes (minimum count of 10 in at least half of samples)

### 2. Metadata Processing
- Parse GEO series matrix file to extract sample metadata
- Standardize group labels: `T1D` or `Healthy`
- Align sample IDs between count matrix and metadata

### 3. Differential Expression Analysis (DESeq2)
- Build a `DESeqDataSet` object using raw counts and sample metadata
- Run DESeq2 normalization and statistical testing (T1D vs. Healthy)
- Filter results: adjusted p-value < 0.05 and |log2FoldChange| > 1

### 4. Visualizations

| Plot | Purpose |
|---|---|
| **MA Plot** | Log2 fold change vs. mean expression across all genes |
| **Enhanced MA Plot** | MA plot with significant genes highlighted in red |
| **Volcano Plot** | Significance vs. magnitude of expression change, colored by direction |
| **Heatmap** | Z-score normalized expression of all 76 DEGs across 46 samples |
| **Bar Plot (DEGs)** | Count of upregulated vs. downregulated genes |
| **Top 10 Genes Bar** | Log2 fold change of the 10 most significant DEGs |
| **PCA Plot** | Sample clustering using variance-stabilized counts (VST) |

### 5. Pathway Enrichment Analysis
- Map significant DEG symbols to Entrez IDs via `bitr()`
- **KEGG enrichment** — run `enrichKEGG()` with adjusted p-value cutoff of 0.05
- **GO enrichment** (Biological Processes) — run `enrichGO()` with adjusted p-value cutoff of 0.05
- Visualize both as dot plots

---

## Results Summary

- **76 DEGs** identified between T1D and healthy whole blood samples
  - **46 upregulated**
  - **30 downregulated**

- **4 proposed biomarker candidates:**

| Gene | Direction | Notes |
|---|---|---|
| `MTRNR2L10` | Upregulated | Novel candidate; limited prior T1D literature |
| `LINC00570` | Upregulated | Novel candidate; limited prior T1D literature |
| `LINC01359` | Downregulated | Previously associated with Hepatocellular Carcinoma |
| `FDPSP2` | Downregulated | Novel candidate; very limited literature overall |

- **Enriched Pathways (GO analysis — 10 pathways identified):**
  - Most enriched: **Cytoplasmic Translation**
  - Also enriched: ribonucleoprotein complex biogenesis, ribosome biogenesis, ribosomal small subunit biogenesis, hydrogen peroxide metabolic process, ribosome assembly, gas transport, hydrogen peroxide catabolic process, oxygen transport, porphyrin-containing compound biosynthetic process

- **KEGG analysis** identified 2 pathways: Ribosome and Coronavirus disease (COVID-19); the GO analysis is considered a more complete representation of enriched biological processes

---

## How to Run

1. Clone this repository:
   ```bash
   git clone https://github.com/raghadmalkawi/Diabetes-type-1-project.git
   cd Diabetes-type-1-project
   ```

2. Download the required data files from [GEO accession GSE123658](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658) and place them in your working directory.

3. Update the file paths at the top of `Type 1 diabetes project.R` to match your local file locations:
   ```r
   raw_data_1 <- read.table(gzfile("path/to/GSE123658_raw_counts_GRCh38.p13_NCBI.tsv.gz"), ...)
   file_path <- "path/to/Human.GRCh38.p13.annot.tsv.gz"
   lines <- readLines(gzfile("path/to/GSE123658-GPL18573_series_matrix.txt.gz"))
   ```

4. Open the script in RStudio and run section by section, or run the entire script at once.

5. Output CSV files and plots will be saved to your working directory (use `getwd()` to confirm the location).
