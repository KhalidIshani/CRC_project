
# ==============================================================================
# Full R Pipeline Rewrite for TCGA-COAD Data Processing,
# CIBERSORT Preparation, and CMS Subtyping
#
# This script incorporates fixes for:
# - Correctly handling gene identifiers (Ensembl to Entrez mapping).
# - Resolving "duplicate row.names" errors.
# - Ensuring proper TPM calculation and log-normalization for CMSclassifier.
# ==============================================================================

# 0. Install and Load Libraries
# It's good practice to put installations at the very top,
# but commented out so they don't run every time.
# Ensure BiocManager is installed first for Bioconductor packages.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install CRAN packages
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) { # DESeq2 is Bioconductor
  BiocManager::install("DESeq2")
}
if (!requireNamespace("immunedeconv", quietly = TRUE)) { # immunedeconv is Bioconductor/GitHub
  BiocManager::install("immunedeconv")
}

# Install Bioconductor packages
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("org.Hs.eg.db") # For human gene annotations

# Install CMSclassifier from GitHub
remotes::install_github("Sage-Bionetworks/CMSclassifier")

# Load all necessary libraries
library(tidyverse) # Includes dplyr, readr, etc.
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2) # Used for normalization in some contexts, but not directly in this TPM calc.
library(immunedeconv) # For CIBERSORT (or other deconvolution methods)
library(org.Hs.eg.db) # For gene ID mapping
library(CMSclassifier) # For human CMS classification
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(scales)
library(ggalluvial)
library(tidyverse)


# --- Data Download and Initial SummarizedExperiment Preparation ---
message("1. Downloading TCGA-COAD data...")
COAD_query <- GDCquery(project = 'TCGA-COAD',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       workflow.type = 'STAR - Counts')

GDCdownload(query = COAD_query)

COAD_data_obj <- GDCprepare(query = COAD_query, summarizedExperiment = TRUE)

# --- Gene Annotations ---
message("2. Extracting gene annotations...")
COAD_genes_meta <- COAD_data_obj@rowRanges %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type, width) # Keep gene_id (Ensembl) as primary key

COAD_genes_meta$width_kb <- COAD_genes_meta$width / 1000

# ================================================================
# --- Robust TPM Conversion and Log2 Transformation ---
# This section is crucial. We will keep Ensembl IDs as row names
# through the TPM calculation, and then map to Entrez IDs for CMS.

message("3. Processing raw counts and converting to log2(TPM + 1)...")
# Start with the assay data, where row names are Ensembl IDs
raw_counts_df <- as.data.frame(assay(COAD_data_obj))

# Add Ensembl IDs as a column for joining with meta data
raw_counts_df <- raw_counts_df %>%
  rownames_to_column("ensembl_id")

# Join with metadata to get gene length (width_kb) and other info.
# We will use 'ensembl_id' as the key throughout.
counts_with_meta <- raw_counts_df %>%
  left_join(COAD_genes_meta %>% dplyr::select(gene_id, width_kb, gene_name, gene_type),
            by = c("ensembl_id" = "gene_id")) %>%
  filter(!is.na(width_kb) & width_kb > 0) # Ensure valid width for TPM calculation

# Identify sample columns (all numeric columns except width_kb)
# Need to exclude gene_name, gene_type, ensembl_id
sample_cols <- names(counts_with_meta)[sapply(counts_with_meta, is.numeric) & !names(counts_with_meta) %in% c("width_kb")]

# Perform TPM calculation using dplyr's across
# This avoids the apply issues and keeps the data frame structure
COAD_tpm_df <- counts_with_meta %>%
  # Divide counts by gene length (kb)
  mutate(across(all_of(sample_cols), ~ . / width_kb)) %>%
  # Normalize to per million
  mutate(across(all_of(sample_cols), ~ (. / sum(., na.rm = TRUE)) * 1e6)) %>%
  # Select only Ensembl ID and TPM values, filter out non-essential info for expression matrix
  dplyr::select(ensembl_id, all_of(sample_cols))

# Set Ensembl IDs as row names for the TPM matrix
COAD_counts_tpm <- COAD_tpm_df %>%
  column_to_rownames("ensembl_id") %>%
  as.matrix() # Convert to matrix as some downstream functions prefer it

# Apply log2 transformation with a pseudocount
COAD_counts_tpm_log2 <- log2(COAD_counts_tpm + 1)

# ================================================================
# --- Sample Filtering and Renaming ---
message("4. Filtering and renaming samples...")

# Load sample lists
# Adjust these paths to your actual file locations
A <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/A.csv")
AK <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AK.csv")
AKP <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKP.csv")
AKPS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKPS.csv")

A_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/A.csv")
AK_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AK_genomic.csv")
AKP_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AKP_genomic.csv")
AKPS_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AKPS_genomic.csv")



A_samples <- A_TCGA$Sample.ID
AK_samples <- AK_TCGA$Sample.ID
AKP_samples <- AKP_TCGA$Sample.ID
AKPS_samples <- AKPS_TCGA$Sample.ID


# Process column names to keep only one sample per patient (based on lexicographical plate ID)
col_info <- data.frame(
  full_id = colnames(COAD_counts_tpm_log2), # Use the processed TPM matrix here
  stringsAsFactors = FALSE
)

col_info$base_id <- substr(col_info$full_id, 1, 15)
col_info$plate <- sub(".*-01R-([^-]+)-.*", "\\1", col_info$full_id)

col_info_unique <- col_info %>%
  group_by(base_id) %>%
  slice_max(order_by = plate, n = 1, with_ties = FALSE) %>%
  ungroup()

# Subset the log2(TPM+1) data to only the selected unique columns
COAD_counts_tpm_log2_unique_samples <- COAD_counts_tpm_log2[, col_info_unique$full_id]

# Rename columns to their 15-character base IDs
colnames(COAD_counts_tpm_log2_unique_samples) <- substr(colnames(COAD_counts_tpm_log2_unique_samples), 1, 15)

# Separate samples into AK and AKP groups
COAD_counts_tpm_AKP <- COAD_counts_tpm_log2_unique_samples[, colnames(COAD_counts_tpm_log2_unique_samples) %in% AKP_samples]
COAD_counts_tpm_AK <- COAD_counts_tpm_log2_unique_samples[, colnames(COAD_counts_tpm_log2_unique_samples) %in% AK_samples]

# Check if any samples were actually selected
if (ncol(COAD_counts_tpm_AKP) == 0) warning("No AKP samples found in the expression data.")
if (ncol(COAD_counts_tpm_AK) == 0) warning("No AK samples found in the expression data.")

# ================================================================
# --- CIBERSORT Preparation ---
# CIBERSORT requires non-log-transformed TPM values with HGNC symbols.
# COAD_counts_tpm_unique_samples has Ensembl IDs as row names and is non-log-transformed TPM.

message("5. Preparing data for CIBERSORT.")

# 1. Get Ensembl IDs for mapping to HGNC
COAD_counts_tpm_unique_samples <- COAD_counts_tpm[, col_info_unique$full_id]

colnames(COAD_counts_tpm_unique_samples) <- substr(colnames(COAD_counts_tpm_unique_samples), 1, 15)

ensembl_ids_for_hgnc <- rownames(COAD_counts_tpm_unique_samples)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids_for_hgnc) # Removes ".XX" suffix

# 2. Map Ensembl IDs to HGNC symbols
message("Mapping Ensembl IDs to HGNC symbols for CIBERSORT input...")
hgnc_symbols_for_cibersort <- mapIds(org.Hs.eg.db,
                                     keys = ensembl_ids_clean,
                                     column = "SYMBOL",     # Request HGNC symbols
                                     keytype = "ENSEMBL",   # Input keys are Ensembl IDs
                                     multiVals = "first")   # Take the first HGNC symbol if multiple mappings exist

# 3. Handle unmapped genes (those that returned NA)
mapped_hgnc_indices <- !is.na(hgnc_symbols_for_cibersort)
if (sum(!mapped_hgnc_indices) > 0) {
  warning(paste0(sum(!mapped_hgnc_indices), " Ensembl IDs could not be mapped to HGNC symbols for CIBERSORT and will be removed."))
  # Optionally, print a sample of unmapped IDs:
  # print(head(ensembl_ids_for_hgnc[!mapped_hgnc_indices], 20))
}

# Filter the TPM data to keep only successfully mapped genes
COAD_tpm_for_cibersort_mapped <- COAD_counts_tpm_unique_samples[mapped_hgnc_indices, ]
hgnc_symbols_mapped <- hgnc_symbols_for_cibersort[mapped_hgnc_indices]

# 4. Handle potential duplicate HGNC symbols by averaging expression values
if (any(duplicated(hgnc_symbols_mapped))) {
  message("Handling duplicate HGNC symbols for CIBERSORT input by averaging expression values...")
  
  df_to_aggregate_hgnc <- COAD_tpm_for_cibersort_mapped %>%
    as.data.frame() %>%
    mutate(HGNC_Symbol = hgnc_symbols_mapped) # Add HGNC Symbol as a column
  
  COAD_tpm_cibersort_final <- df_to_aggregate_hgnc %>%
    group_by(HGNC_Symbol) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    as.data.frame()
  
  rownames(COAD_tpm_cibersort_final) <- COAD_tpm_cibersort_final$HGNC_Symbol
  COAD_tpm_cibersort_final$HGNC_Symbol <- NULL
  
} else {
  COAD_tpm_cibersort_final <- COAD_tpm_for_cibersort_mapped
  rownames(COAD_tpm_cibersort_final) <- hgnc_symbols_mapped
}

message("CIBERSORT input data prepared with HGNC symbols as row names.")

# Now, subset this HGNC-symbol-mapped TPM data into AK and AKP groups
COAD_counts_tpm_AK_for_cibersort <- COAD_tpm_cibersort_final[,colnames(COAD_tpm_cibersort_final) %in% AK_samples]
COAD_counts_tpm_AKP_for_cibersort <- COAD_tpm_cibersort_final[,colnames(COAD_tpm_cibersort_final) %in% AKP_samples]
COAD_counts_tpm_A_for_cibersort <- COAD_tpm_cibersort_final[,colnames(COAD_tpm_cibersort_final) %in% A_samples]
COAD_counts_tpm_AKPS_for_cibersort <- COAD_tpm_cibersort_final[,colnames(COAD_tpm_cibersort_final) %in% AKPS_samples]

# Add GeneSymbol column for CIBERSORT, which typically expects it
COAD_counts_tpm_AK_for_cibersort <- COAD_counts_tpm_AK_for_cibersort %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneSymbol")

COAD_counts_tpm_AKP_for_cibersort <- COAD_counts_tpm_AKP_for_cibersort %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneSymbol")

COAD_counts_tpm_A_for_cibersort <- COAD_counts_tpm_A_for_cibersort %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneSymbol")

COAD_counts_tpm_AKPS_for_cibersort <- COAD_counts_tpm_AKPS_for_cibersort %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneSymbol")

# Write to file for CIBERSORT. Adjust paths as necessary.
write.table(COAD_counts_tpm_AK_for_cibersort,
            file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/COAD_counts_tpm_AK.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(COAD_counts_tpm_AKP_for_cibersort,
            file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/COAD_counts_tpm_AKP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(COAD_counts_tpm_A_for_cibersort,
            file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/COAD_counts_tpm_A.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(COAD_counts_tpm_AKPS_for_cibersort,
            file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/COAD_counts_tpm_AKPS.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
# ================================================================
# --- Statistical Analysis of CIBERSORT Results ---
message("6. Performing statistical analysis on CIBERSORT results- A vs. AK (assuming external run).")

CIBERSORT_A_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_A_immunefraction.csv")
CIBERSORT_AK_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AK_immunefraction.csv")
results_A_vs_AK <- lapply(colnames(CIBERSORT_A_RESULTS[,2:23]), function(cell_type) {
  group1 <- CIBERSORT_A_RESULTS[[cell_type]] # Data for group A
  group2 <- CIBERSORT_AK_RESULTS[[cell_type]] # Data for group AK
  
  test <- tryCatch(wilcox.test(group1, group2, paired = FALSE, alternative = "less"),
                   error = function(e) {
                     warning(paste("Wilcoxon test failed for cell type", cell_type, " (A vs AK):", e$message))
                     return(list(p.value = NA)) # Return NA if test fails
                   })
  
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  
  return(c(
    mean_A = mean1,    # Renamed for A group
    mean_AK = mean2,   # Renamed for AK group
    p_value = test$p.value
  ))
})

results_df_A_vs_AK <- as.data.frame(do.call(rbind, results_A_vs_AK))

# Add cell type names and reorder columns
results_df_A_vs_AK$cell_type <- colnames(CIBERSORT_A_RESULTS[,2:23])
results_df_A_vs_AK <- results_df_A_vs_AK[,c("cell_type", "mean_A", "mean_AK", "p_value")]

# Apply FDR correction (Benjamini-Hochberg)
results_df_A_vs_AK$FDR_corrected_p_value <- p.adjust(results_df_A_vs_AK$p_value, method="BH")

results_df_A_vs_AK <- results_df_A_vs_AK[order(results_df_A_vs_AK$p_value), ]

# Write statistical results for A vs AK comparison
write.csv(results_df_A_vs_AK,
          "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_results_A_vs_AK.csv",
          row.names = FALSE)

message("CIBERSORT A vs AK comparison complete.")


message("6. Performing statistical analysis on CIBERSORT results- AK vs. AKP (assuming external run).")
# Load CIBERSORT results. Adjust paths.
CIBERSORT_AK_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AK_immunefraction.csv")
CIBERSORT_AKP_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKP_immunefraction.csv")

# Perform Wilcoxon rank-sum test for each cell type
results <- lapply(colnames(CIBERSORT_AK_RESULTS[,2:23]), function(cell_type) {
  group1 <- CIBERSORT_AK_RESULTS[[cell_type]]
  group2 <- CIBERSORT_AKP_RESULTS[[cell_type]]
  
  test <- tryCatch(wilcox.test(group1, group2, paired = FALSE, alternative = "less"),
                   error = function(e) {
                     warning(paste("Wilcoxon test failed for cell type", cell_type, ":", e$message))
                     return(list(p.value = NA)) # Return NA if test fails
                   })
  
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  
  return(c(
    mean_AK = mean1,
    mean_AKP = mean2,
    p_value = test$p.value
  ))
})

results_df <- as.data.frame(do.call(rbind, results))

# Add cell type names and reorder columns
results_df$cell_type <- colnames(CIBERSORT_AK_RESULTS[,2:23])
results_df <- results_df[,c("cell_type", "mean_AK", "mean_AKP", "p_value")]

# Apply FDR correction (Benjamini-Hochberg)
results_df$FDR_corrected_p_value <- p.adjust(results_df$p_value, method="BH")

results_df <- results_df[order(results_df$p_value), ]

# Write statistical results
write.csv(results_df, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_results_AK_vs_AKP.csv", row.names = FALSE)
message("CIBERSORT AK vs AKP comparison complete.")


message("Performing statistical analysis: CIBERSORT_AKP_RESULTS vs. CIBERSORT_AKPS_RESULTS")

CIBERSORT_AKP_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKP_immunefraction.csv")
CIBERSORT_AKPS_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKPS_immunefraction.csv")
# Assuming CIBERSORT_AKPS_RESULTS is loaded from a file like:
# CIBERSORT_AKPS_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKPS_immunefraction.csv") # Adjust path

results_AKP_vs_AKPS <- lapply(colnames(CIBERSORT_AKP_RESULTS[,2:23]), function(cell_type) { # Use AKP for column names reference
  group1 <- CIBERSORT_AKP_RESULTS[[cell_type]]  # Data for group AKP
  group2 <- CIBERSORT_AKPS_RESULTS[[cell_type]] # Data for group AKPS
  
  test <- tryCatch(wilcox.test(group1, group2, paired = FALSE, alternative = "less"),
                   error = function(e) {
                     warning(paste("Wilcoxon test failed for cell type", cell_type, " (AKP vs AKPS):", e$message))
                     return(list(p.value = NA)) # Return NA if test fails
                   })
  
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  
  return(c(
    mean_AKP = mean1,    # Renamed for AKP group
    mean_AKPS = mean2,   # Renamed for AKPS group
    p_value = test$p.value
  ))
})

results_df_AKP_vs_AKPS <- as.data.frame(do.call(rbind, results_AKP_vs_AKPS))

# Add cell type names and reorder columns
results_df_AKP_vs_AKPS$cell_type <- colnames(CIBERSORT_AKP_RESULTS[,2:23]) # Use AKP for column names reference
results_df_AKP_vs_AKPS <- results_df_AKP_vs_AKPS[,c("cell_type", "mean_AKP", "mean_AKPS", "p_value")]

# Apply FDR correction (Benjamini-Hochberg)
results_df_AKP_vs_AKPS$FDR_corrected_p_value <- p.adjust(results_df_AKP_vs_AKPS$p_value, method="BH")

results_df_AKP_vs_AKPS <- results_df_AKP_vs_AKPS[order(results_df_AKP_vs_AKPS$p_value), ]

# Write statistical results for AKP vs AKPS comparison
write.csv(results_df_AKP_vs_AKPS,
          "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_results_AKP_vs_AKPS.csv",
          row.names = FALSE)

message("CIBERSORT AKP vs AKPS comparison complete.")





####make alluvial plot of the immune subtype data ######

read_cibersort <- function(file, group_name) {
  df <- read.csv(file)
  
  # Keep only immune fractions (exclude Mixture + QC cols)
  df <- df %>%
    dplyr::select(-Mixture, -P.value, -Correlation, -RMSE, -`Absolute.score..sig.score.`)
  
  # Collapse subtypes into parent categories
  df_collapsed <- df %>%
    transmute(
      B_cells = `B.cells.naive` + `B.cells.memory`,
      Plasma_cells = `Plasma.cells`,
      T_cells = `T.cells.CD8` + `T.cells.CD4.naive` +
        `T.cells.CD4.memory.resting` + `T.cells.CD4.memory.activated` +
        `T.cells.follicular.helper` + `T.cells.regulatory..Tregs.` +
        `T.cells.gamma.delta`,
      NK_cells = `NK.cells.resting` + `NK.cells.activated`,
      Monocytes = `Monocytes`,
      Macrophages = `Macrophages.M0` + `Macrophages.M1` + `Macrophages.M2`,
      Dendritic_cells = `Dendritic.cells.resting` + `Dendritic.cells.activated`,
      Mast_cells = `Mast.cells.resting` + `Mast.cells.activated`,
      Eosinophils = `Eosinophils`,
      Neutrophils = `Neutrophils`
    )
  
  df_collapsed$Group <- group_name
  return(df_collapsed)
}


A <- read_cibersort("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_A_immunefraction.csv", "A")
AK <- read_cibersort("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AK_immunefraction.csv", "AK")
AKP <- read_cibersort("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKP_immunefraction.csv", "AKP")
AKPS <- read_cibersort("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKPS_immunefraction.csv", "AKPS")

df <- bind_rows(A, AK, AKP, AKPS)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Group,
    names_to = "CellType",
    values_to = "Fraction"
  )

# Average per group
df_summary <- df_long %>%
  group_by(Group, CellType) %>%
  summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

# Order stages
df_summary$Group <- factor(df_summary$Group, levels = c("A", "AK", "AKP", "AKPS"))

# Plot alluvial

humanAlluvial <- ggplot(df_summary,
       aes(x = Group, stratum = CellType, alluvium = CellType,
           y = Fraction, fill = CellType)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.7) +
  geom_stratum() +
  theme_minimal() +
  labs(title = "Immune Cell Fractions Across Tumor Progression",
       y = "Fraction of Immune Cells", x = "Progression Stage") +
  theme(legend.position = "right")


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/Human_alluvialplot.pdf",humanAlluvial, width = 14, height = 8)

##########################
# ================================================================
# --- CMS Subtyping for Human Samples ---
message("7. Performing CMS subtyping for human samples.")

# The input for CMSclassifier needs to be log2(TPM+1) with Entrez IDs as row names.
# We will use COAD_counts_tpm_log2_unique_samples which has Ensembl IDs as row names
# before subsetting into AK and AKP. This ensures we process all relevant genes.

# 1. Get Ensembl IDs from the processed log2(TPM+1) data
ensembl_ids_for_cms <- rownames(COAD_counts_tpm_log2_unique_samples)
ensembl_ids_clean_cms <- sub("\\..*$", "", ensembl_ids_for_cms) # Remove version numbers

# 2. Map Ensembl IDs to Entrez IDs
message("Mapping Ensembl IDs to Entrez IDs for CMS classification...")
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids_clean_cms,
                     column = "ENTREZID",
                     keytype = "ENSEMBL", # CRITICAL: Use ENSEMBL as keytype
                     multiVals = "first") # Take the first Entrez ID if multiple mappings exist

# 3. Handle unmapped genes (those that returned NA)
mapped_indices <- !is.na(entrez_ids)
if (sum(!mapped_indices) > 0) {
  warning(paste0(sum(!mapped_indices), " Ensembl IDs could not be mapped to Entrez IDs for CMS and will be removed."))
  # Optionally, print a sample of unmapped IDs for inspection:
  # print(head(ensembl_ids_for_cms[!mapped_indices], 20))
}

# Filter the expression data to keep only successfully mapped genes
COAD_CMS_input_mapped <- COAD_counts_tpm_log2_unique_samples[mapped_indices, ]
entrez_ids_mapped <- entrez_ids[mapped_indices]

# 4. Handle potential duplicate Entrez IDs by averaging expression values
if (any(duplicated(entrez_ids_mapped))) {
  message("Handling duplicate Entrez IDs for CMS input by averaging expression values...")
  
  df_to_aggregate_cms <- COAD_CMS_input_mapped %>%
    as.data.frame() %>%
    mutate(EntrezID = entrez_ids_mapped)
  
  COAD_CMS_input_final <- df_to_aggregate_cms %>%
    group_by(EntrezID) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    as.data.frame()
  
  rownames(COAD_CMS_input_final) <- COAD_CMS_input_final$EntrezID
  COAD_CMS_input_final$EntrezID <- NULL
  
} else {
  COAD_CMS_input_final <- COAD_CMS_input_mapped
  rownames(COAD_CMS_input_final) <- entrez_ids_mapped
}

message("CMS input data prepared with Entrez IDs as row names.")

# Now, run CMS classification
# Note: CMSclassifier expects samples in columns and Entrez IDs in rows.
# COAD_CMS_input_final is in this format.

message("Running CMS classification using CMSclassifier...")
cms_results <- classifyCMS(as.matrix(COAD_CMS_input_final), method = "RF") # 'method = "RF"' is robust

# Print CMS classification results
message("CMS Classification Results:")
print(cms_results$predictedCMS) # Or whatever the primary results object is named in your CMSclassifier version

# You can also save these results to a file
write.csv(cms_results$predictedCMS,
          "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classification_results.csv",
          row.names = TRUE)

message("Pipeline execution complete.")


CMS_classifications <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classification_results.csv")


CMS_classiciations_AK <- CMS_classifications[rownames(CMS_classifications) %in% AK_samples, , drop=FALSE]
CMS_classiciations_AKP <- CMS_classifications[rownames(CMS_classifications) %in% AKP_samples, , drop=FALSE]
CMS_classiciations_A <- CMS_classifications[rownames(CMS_classifications) %in% A_samples, , drop=FALSE]
CMS_classiciations_AKPS <- CMS_classifications[rownames(CMS_classifications) %in% AKPS_samples, , drop=FALSE]
#### now try running CIBERSORT on the individual CMS subtypes comparing them ####### 

message("Starting CIBERSORT analysis comparing multiple groups within each CMS subtype.")

# Ensure CMS_classifications is correctly initialized from cms_results$predictedCMS
# and has 'Sample Name' as a column (as per previous fixes).
# Example:
 CMS_classifications <- cms_results$predictedCMS %>%
   rownames_to_column("Sample Name") %>%
   as.data.frame()

 
 
# Load all CIBERSORT results dataframes (adjust paths as needed)
# CIBERSORT_A_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_A_immunefraction.csv")
# CIBERSORT_AK_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AK_immunefraction.csv")
# CIBERSORT_AKP_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKP_immunefraction.csv")
# CIBERSORT_AKPS_RESULTS <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_AKPS_immunefraction.csv")


# Define the comparisons you want to make
# Each element in this list is a vector of two strings: c("Group1_Results_DF_Name", "Group2_Results_DF_Name")
# And their corresponding mean column names for the output dataframe
message("Starting CIBERSORT analysis comparing multiple groups WITHIN each CMS subtype.")

# Loop through each CMS subtype and perform the analysis
cms_subtypes <- c("CMS1", "CMS2", "CMS3", "CMS4")

# Define the comparisons you want to make within each CMS subtype
comparison_pairs_within_cms <- list(
  list(df1_name = "CIBERSORT_A_RESULTS", df2_name = "CIBERSORT_AK_RESULTS",
       mean1_col = "mean_A", mean2_col = "mean_AK", filename_suffix = "A_vs_AK"),
  list(df1_name = "CIBERSORT_AK_RESULTS", df2_name = "CIBERSORT_AKP_RESULTS",
       mean1_col = "mean_AK", mean2_col = "mean_AKP", filename_suffix = "AK_vs_AKP"),
  list(df1_name = "CIBERSORT_AKP_RESULTS", df2_name = "CIBERSORT_AKPS_RESULTS",
       mean1_col = "mean_AKP", mean2_col = "mean_AKPS", filename_suffix = "AKP_vs_AKPS")
)

for (cms_type in cms_subtypes) {
  message(paste0("\nProcessing CMS type: ", cms_type))
  
  # Filter samples for the current CMS type
  current_cms_samples <- CMS_classifications %>%
    dplyr::filter(RF == cms_type) # Use the CMS_classifications with 'Sample Name' column
  
  # Check if there are enough samples in this CMS subtype for any analysis
  if (nrow(current_cms_samples) == 0) {
    warning(paste("No samples found for CMS type", cms_type, ". Skipping all comparisons for this subtype."))
    next # Skip to the next CMS type in the loop
  }
  
  # Iterate through each defined comparison pair
  for (comp_pair in comparison_pairs_within_cms) {
    df1_name <- comp_pair$df1_name
    df2_name <- comp_pair$df2_name
    mean1_col <- comp_pair$mean1_col
    mean2_col <- comp_pair$mean2_col
    filename_suffix <- comp_pair$filename_suffix
    
    message(paste0("  Comparing ", df1_name, " vs ", df2_name, " within CMS ", cms_type))
    
    # Get the actual dataframes using get()
    df1_results_full <- get(df1_name)
    df2_results_full <- get(df2_name)
    
    # Subset CIBERSORT results for the current CMS type and desired comparison groups
    df1_results_cms_subset <- df1_results_full %>%
      dplyr::filter(Mixture %in% current_cms_samples$`Sample Name`)
    
    df2_results_cms_subset <- df2_results_full %>%
      dplyr::filter(Mixture %in% current_cms_samples$`Sample Name`)
    
    # --- Diagnostic lines for "all NA values" issue ---
    message(paste0("  CMS Classification Samples (current CMS type): ", nrow(current_cms_samples)))
    # print(head(current_cms_samples$`Sample Name`, 5)) # Uncomment to inspect
    
    message(paste0("  Full CIBERSORT ", df1_name, " samples: ", nrow(df1_results_full)))
    # print(head(df1_results_full$Mixture, 5)) # Uncomment to inspect
    
    message(paste0("  Full CIBERSORT ", df2_name, " samples: ", nrow(df2_results_full)))
    # print(head(df2_results_full$Mixture, 5)) # Uncomment to inspect
    
    overlap1 <- sum(df1_results_full$Mixture %in% current_cms_samples$`Sample Name`)
    overlap2 <- sum(df2_results_full$Mixture %in% current_cms_samples$`Sample Name`)
    message(paste0("  Overlap between ", df1_name, " Mixture and CMS samples: ", overlap1, " (expected >0)"))
    message(paste0("  Overlap between ", df2_name, " Mixture and CMS samples: ", overlap2, " (expected >0)"))
    
    message(paste0("  After filtering - ", df1_name, " subsetted samples: ", nrow(df1_results_cms_subset)))
    message(paste0("  After filtering - ", df2_name, " subsetted samples: ", nrow(df2_results_cms_subset)))
    # --- End Diagnostic lines ---
    
    # Check for sufficient samples in each group for Wilcoxon test
    if (nrow(df1_results_cms_subset) < 2 || nrow(df2_results_cms_subset) < 2) {
      warning(paste0("  Not enough samples (at least 2 per group) for Wilcoxon test in CMS ", cms_type, " (", filename_suffix, "). Skipping statistical analysis."))
      
      # Prepare an empty results_df for consistency
      results_df <- data.frame(cell_type = colnames(df1_results_full[,2:23]),
                               p_value = NA, FDR_corrected_p_value = NA)
      results_df[[mean1_col]] <- NA
      results_df[[mean2_col]] <- NA
      results_df <- results_df[, c("cell_type", mean1_col, mean2_col, "p_value", "FDR_corrected_p_value")]
      
    } else {
      # Perform Wilcoxon test for each cell type
      results_list <- lapply(colnames(df1_results_cms_subset[,2:23]), function(cell_type) {
        group1 <- df1_results_cms_subset[[cell_type]]
        group2 <- df2_results_cms_subset[[cell_type]]
        
        # Ensure there's data in both groups for the test, otherwise wilcox.test fails
        if (length(group1[!is.na(group1)]) == 0 || length(group2[!is.na(group2)]) == 0) {
          p_val <- NA
          warning(paste("  Not enough non-NA values for Wilcoxon test for cell type", cell_type, "in CMS", cms_type, " (", filename_suffix, ")"))
        } else {
          test <- tryCatch(wilcox.test(group1, group2, paired = FALSE, alternative = "less"),
                           error = function(e) {
                             warning(paste("  Wilcoxon test failed for cell type", cell_type, "in CMS", cms_type, " (", filename_suffix, "):", e$message))
                             return(list(p.value = NA))
                           })
          p_val <- test$p.value
        }
        
        mean1 <- mean(group1, na.rm = TRUE)
        mean2 <- mean(group2, na.rm = TRUE)
        
        # Return a named vector with dynamic mean column names
        res_vec <- c(p_value = p_val)
        res_vec[mean1_col] <- mean1
        res_vec[mean2_col] <- mean2
        return(res_vec)
      })
      
      results_df <- as.data.frame(do.call(rbind, results_list))
      
      # Add cell type names and reorder columns dynamically
      results_df$cell_type <- colnames(df1_results_cms_subset[,2:23])
      results_df <- results_df[,c("cell_type", mean1_col, mean2_col, "p_value")]
      
      # Apply FDR correction (Benjamini-Hochberg)
      results_df$FDR_corrected_p_value <- p.adjust(results_df$p_value, method="BH")
      
      results_df <- results_df[order(results_df$p_value), ]
    }
    
    # Write statistical results for the current CMS type and comparison
    write.csv(results_df,
              file = paste0("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/CIBERSORT_results_CMS", cms_type, "_", filename_suffix, ".csv"),
              row.names = FALSE)
  }
}

message("CIBERSORT analysis by CMS subtype for all specified comparisons complete.")

##### now that there is no statistical significance found, let's try GSEA to see if there are pathway differences 
# 1. Start with raw counts matrix from COAD_data_obj
raw_counts <- assay(COAD_data_obj)  # genes x full TCGA sample names

# 2. Build sample info table
col_info <- data.frame(
  full_id = colnames(raw_counts),
  stringsAsFactors = FALSE
)

# 3. Extract base ID (first 15 characters of TCGA barcode)
col_info$base_id <- substr(col_info$full_id, 1, 15)

# 4. Extract plate ID (helps select the latest replicate)
col_info$plate <- sub(".*-01R-([^-]+)-.*", "\\1", col_info$full_id)

# 5. Keep only the latest sample per base_id (if multiple replicates exist)
library(dplyr)

col_info_unique <- col_info %>%
  group_by(base_id) %>%
  slice_max(order_by = plate, n = 1, with_ties = FALSE) %>%
  ungroup()

# 6. Subset raw counts to only these selected unique samples
raw_counts_unique <- raw_counts[, col_info_unique$full_id]

# 7. Rename columns to 15-character base IDs (matches your Excel lists)
colnames(raw_counts_unique) <- substr(colnames(raw_counts_unique), 1, 15)

# Label sample metadata
sample_ids <- colnames(raw_counts_unique)

sample_metadata <- data.frame(
  sample_id = sample_ids,
  group = case_when(
    sample_ids %in% A_samples ~ "A",
    sample_ids %in% AK_samples ~ "AK",
    sample_ids %in% AKP_samples ~ "AKP",
    sample_ids %in% AKPS_samples ~ "AKPS",
    TRUE ~ NA_character_
  )
) %>% filter(!is.na(group))

dds <- DESeqDataSetFromMatrix(countData = raw_counts_unique[, sample_metadata$sample_id],
                              colData = sample_metadata,
                              design = ~ group)

dds <- DESeq(dds)

comparisons <- list(
  c("AK", "A"),
  c("AKP", "AK"),
  c("AKPS", "AKP"))

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA")

for (cmp in comparisons) {
  res <- results(dds, contrast = c("group", cmp[1], cmp[2]))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(COAD_genes_meta %>% dplyr::select(gene_id, gene_name), 
              by = c("ensembl_id" = "gene_id")) %>%
    filter(!is.na(stat) & !is.na(gene_name)) %>%
    arrange(desc(stat)) %>%
    dplyr::select(gene_name, stat)
  
  fname <- paste0(cmp[1], "_vs_", cmp[2], "_ranked_genes.rnk")
  write.table(res_df, file = fname, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}
########### now that GSEA has been run -- make a heat map ######### 
AvAK_neg <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AvsAK.GseaPreranked.1756686396995/gsea_report_for_na_neg_1756686396995.tsv")
AvAK_pos <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AvsAK.GseaPreranked.1756686396995/gsea_report_for_na_pos_1756686396995.tsv")
AvAK <- rbind(AvAK_neg,AvAK_pos)
names(AvAK)[names(AvAK) == "NAME"] <- "pathway"
names(AvAK)[names(AvAK) == "FDR q-val"] <- "padj"
AvAK$NES <- as.numeric(AvAK$NES)

AKvAKP_neg <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AKvsAKP.GseaPreranked.1756686558945/gsea_report_for_na_neg_1756686558945.tsv")
AKvAKP_pos <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AKvsAKP.GseaPreranked.1756686558945/gsea_report_for_na_pos_1756686558945.tsv")
AKvAKP <- rbind(AKvAKP_neg,AKvAKP_pos)
names(AKvAKP)[names(AKvAKP) == "NAME"] <- "pathway"
names(AKvAKP)[names(AKvAKP) == "FDR q-val"] <- "padj"
AKvAKP$NES <- as.numeric(AKvAKP$NES)

AKPvAKPS_neg <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AKPvsAKPS.GseaPreranked.1756686601595/gsea_report_for_na_neg_1756686601595.tsv")
AKPvAKPS_pos <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Output/AKPvsAKPS.GseaPreranked.1756686601595/gsea_report_for_na_pos_1756686601595.tsv")
AKPvAKPS <- rbind(AKPvAKPS_neg,AKPvAKPS_pos)
names(AKPvAKPS)[names(AKPvAKPS) == "NAME"] <- "pathway"
names(AKPvAKPS)[names(AKPvAKPS) == "FDR q-val"] <- "padj"
AKPvAKPS$NES <- as.numeric(AKPvAKPS$NES)

#####prepare heatmap ######

plot_full_hallmark_heatmap_with_priority <- function(gsea_results_list, comparisons, immune_pathways, pval_cutoff = 0.05) {
  # gsea_results_list: named list of fgsea outputs (each a data.frame)
  # comparisons: named vector like c("AK_vs_A" = "A → AK", ...)
  # immune_pathways: character vector of Hallmark pathway names to put at top
  
  # 1. Combine all GSEA results
  all_results <- purrr::imap_dfr(gsea_results_list, ~ {
    .x %>%
      dplyr::select(pathway, NES, padj) %>%
      mutate(comparison = .y)
  })
  
  # 2. Filter for Hallmark pathways
  all_results <- all_results %>%
    filter(grepl("^HALLMARK_", pathway))
  
  # 3. Create NES matrix
  nes_matrix <- all_results %>%
    dplyr::select(pathway, comparison, NES) %>%
    pivot_wider(names_from = comparison, values_from = NES) %>%
    column_to_rownames("pathway")
  
  # 4. Create significance ("+") matrix
  plus_matrix <- all_results %>%
    mutate(sig = ifelse(padj < pval_cutoff, "+", "")) %>%
    dplyr::select(pathway, comparison, sig) %>%
    pivot_wider(names_from = comparison, values_from = sig) %>%
    column_to_rownames("pathway")
  
  # 5. Ensure row order: immune pathways on top, rest alphabetically
  all_pathways <- rownames(nes_matrix)
  immune_order <- immune_pathways[immune_pathways %in% all_pathways]
  remaining <- setdiff(sort(all_pathways), immune_order)
  final_order <- c(immune_order, remaining)
  
  nes_matrix <- nes_matrix[final_order, ]
  plus_matrix <- plus_matrix[final_order, ]
  
  # 6. Clean column names to user-friendly progression names
  colnames(nes_matrix) <- comparisons[colnames(nes_matrix)]
  colnames(plus_matrix) <- comparisons[colnames(plus_matrix)]
  
  # 7. Plot full heatmap
  pheatmap(nes_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",
           display_numbers = plus_matrix,
           number_color = "black",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Hallmark Pathway Enrichment",
           fontsize_row = 9,
           fontsize_col = 9, 
           angle_col = 0)
}


immune_pathways <- c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
  "HALLMARK_TGF_BETA_SIGNALING"
)

oncogenic_pathways <- c(
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)


gsea_results_list <- list(
  "AK_vs_A" = AvAK,
  "AKP_vs_AK" = AKvAKP,
  "AKPS_vs_AKP" = AKPvAKPS
)

# Named vector to control column order and labels
comparisons <- c("AK_vs_A" = "A vs. AK",
                 "AKP_vs_AK" = "AK vs. AKP",
                 "AKPS_vs_AKP"= "AKP vs. AKPS")

gsea_results_list_immmune_AvAK <- gsea_results_list$AK_vs_A[gsea_results_list$AK_vs_A$pathway %in% immune_pathways,]
gsea_results_list_oncogenic_AvAK <- gsea_results_list$AK_vs_A[gsea_results_list$AK_vs_A$pathway %in% oncogenic_pathways,]

gsea_results_list_immmune_AKvAKP <- gsea_results_list$AKP_vs_AK[gsea_results_list$AKP_vs_AK$pathway %in% immune_pathways,]
gsea_results_list_oncogenic_AKvAKP <- gsea_results_list$AKP_vs_AK[gsea_results_list$AKP_vs_AK$pathway %in% oncogenic_pathways,]

gsea_results_list_immmune_AKPvAKPS <- gsea_results_list$AKPS_vs_AKP[gsea_results_list$AKPS_vs_AKP$pathway %in% immune_pathways,]
gsea_results_list_oncogenic_AKPvAKPS <- gsea_results_list$AKPS_vs_AKP[gsea_results_list$AKPS_vs_AKP$pathway %in% oncogenic_pathways,]


gsea_results_list_immune <- list(
  "AK_vs_A"   = gsea_results_list_immmune_AvAK,
  "AKP_vs_AK" = gsea_results_list_immmune_AKvAKP,
  "AKPS_vs_AKP" = gsea_results_list_immmune_AKPvAKPS
)


gsea_results_list_oncogenic <- list(
  "AK_vs_A"   = gsea_results_list_oncogenic_AvAK,
  "AKP_vs_AK" = gsea_results_list_oncogenic_AKvAKP,
  "AKPS_vs_AKP" = gsea_results_list_oncogenic_AKPvAKPS
)


hallMarkHeatMapHuman_IMMUNE  <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list_immune,
  comparisons,
  immune_pathways,
  pval_cutoff = 0.05
)

hallMarkHeatMapHuman_ONCOGENIC  <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list_oncogenic,
  comparisons,
  immune_pathways,
  pval_cutoff = 0.05
)


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/hallMarkHeatMapHuman_immune.pdf",hallMarkHeatMapHuman_IMMUNE, width = 14, height = 8)
ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/hallMarkHeatMapHuman_oncogenic.pdf",hallMarkHeatMapHuman_ONCOGENIC, width = 14, height = 8)



###now make alluvial plot showing different CMS classifications 
CMS_classifications_barplot <- CMS_classifications %>% mutate(group=case_when(
  `Sample Name` %in% A_samples ~ "A", 
  `Sample Name` %in% AK_samples ~ "AK", 
  `Sample Name` %in% AKP_samples ~ "AKP", 
  `Sample Name` %in% AKPS_samples ~ "AKPS", 
))

CMS_classifications_barplot <- CMS_classifications_barplot[!is.na(CMS_classifications_barplot$group),]

CMS_classifications_barplot$group <- factor(CMS_classifications_barplot$group, levels = c("A","AK","AKP","AKPS"))

write.csv(CMS_classifications_barplot, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classifications_alluvial.csv")

CMS_classifications_barplot <- CMS_classifications_barplot[!is.na(CMS_classifications_barplot$RF),]

CMS_classifications_barplot_grouped <- CMS_classifications_barplot %>%
  group_by(RF, group) %>%
  summarise(Freq =n(), .groups = "drop")

CMS_classifications_barplot_grouped <- CMS_classifications_barplot_grouped %>% group_by(group) %>% mutate(
  proportion = Freq/sum(Freq) 
) %>% ungroup()

#make alluvial 

humanAlluvialCMS <- ggplot(CMS_classifications_barplot_grouped,
                        aes(x = group, stratum = RF, alluvium = RF,
                            y = proportion, fill = RF)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.7) +
  geom_stratum() +
  theme_minimal() +
  labs(title = "CMS Subtype Composition Across Tumor Progression",
       y = "Fraction of Samples Belonging to CMS Subtype", x = "Progression Stage") +
  theme(legend.position = "right")


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/Human_alluvialplot_CMS.pdf",humanAlluvialCMS, width = 14, height = 8)

###make PCA of A, AK, AKP, AKPS 

vsd <- vst(dds)
data_human <- as.data.frame(assay(vsd))
pca_result_human <- prcomp(t(data_human))
pca_df_human <- as.data.frame(pca_result_human$x)
pca_df_human <- pca_df_human %>%
  rownames_to_column("sample")
variance_explained <- pca_result_human$sdev^2 / sum(pca_result_human$sdev^2) * 100

pca_df_human <- pca_df_human %>% mutate(grp_gene=case_when(
  sample %in% A_samples ~ "A", 
  sample %in% AK_samples ~ "AK", 
  sample %in% AKP_samples ~ "AKP", 
  sample %in% AKPS_samples ~ "AKPS", 
))


pca_human <- ggplot(pca_df_human, aes(x = PC1, y = PC2, color=grp_gene)) +
  geom_point(size=7) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(variance_explained[1], 2), "% variance)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "% variance)"))  +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),
        axis.ticks = element_line(linewidth = 0.5)) 

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/pca_human.pdf",pca_human , width = 10, height = 8)

#####################################################




#### make umap including both TCGA and CPTAC-2 datasets. Will standardize to z-scores so that the CPTAC-2 data can be integrated. 
library(uwot)
library(ggplot2)
library(dplyr)
library(readr)
library(ComplexUpset)
library(sva)
library(biomaRt)
library(dplyr)
set.seed(123)
library(openxlsx)

CPTAC_2_zscores <- readr::read_delim("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_cptac_2019/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", delim = "\t") 
TCGA_zscores <- readr::read_delim("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt", delim = "\t") 

CPTAC_2_zscores <- CPTAC_2_zscores %>%
  filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Entrez_Gene_Id")

TCGA_zscores <- TCGA_zscores %>%
  filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Entrez_Gene_Id")




colnames(TCGA_zscores) <- sub("([A-Za-z])$", "", colnames(TCGA_zscores)) 

#A_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/A.csv")
AK_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AK.csv")
AKP_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKP.csv")
AKPS_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKPS.csv")

#A_CPTAC <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/A_CPTAC2.csv")
#AK_CPTAC <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AK_CPTAC2.csv")
#AKP_CPTAC <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKP_CPTAC2.csv")




AK_samples <- AK_TCGA$Sample.ID
AKP_samples <- AKP_TCGA$Sample.ID
AKPS_samples <- AKPS_TCGA$Sample.ID
All_defined_samples <- union(AK_samples, AKP_samples)
All_defined_samples <- union(All_defined_samples, AKPS_samples)


No_defined_group <- setdiff(colnames(TCGA_zscores), All_defined_samples)

all_groups <- list(
  AK = AK_samples,
  AKP = AKP_samples,
  AKPS = AKPS_samples,
  Other = No_defined_group
)

sample_to_group <- stack(all_groups) %>%
  dplyr::rename(sample = values, group = ind)


samples_with_group <- sample_to_group$sample


combined_mat_filt  <- TCGA_zscores[, colnames(TCGA_zscores) %in% samples_with_group]

combined_mat_filt <- combined_mat_filt[complete.cases(combined_mat_filt),]


#TCGA <- sum(grepl("^TCGA", colnames(combined_mat_filt)))
#CPTAC <- ncol(combined_mat_filt) - TCGA

#batch <- c(rep("CPTAC", CPTAC), rep("TCGA", TCGA))
#mod <- model.matrix(~ 1, data = data.frame(batch)) 

#batch_corrected_combined_mat_filt <-ComBat(
  #dat = combined_mat_filt,
  #batch = batch,
  #mod = mod,
  #par.prior = TRUE,
  #prior.plots = FALSE
#)

umap_res <- umap(t(combined_mat_filt))

umap_df <- as.data.frame(umap_res$layout)
umap_df$sample <- colnames(combined_mat_filt)

umap_df <- umap_df %>%
  left_join(sample_to_group, by = "sample")


Umap <- ggplot(umap_df, aes(V1, V2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  
  # 95% CI ellipses
  stat_ellipse(aes(fill = group),
               type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.15, color = NA) +
  
  scale_color_manual(
    values = c(
      "AK"   = "#E41A1C",
      "AKP"  = "#377EB8",
      "AKPS" = "#4DAF4A",
      "Other"= "grey70"
    )
  ) +
  scale_fill_manual(
    values = c(
      "AK"   = "#E41A1C",
      "AKP"  = "#377EB8",
      "AKPS" = "#4DAF4A",
      "Other"= "grey70"
    )
  ) +
  
  theme_classic() +
  labs(color = "Mutation Group", fill = "Mutation Group", x="UMAP1", y="UMAP2") +
  
  # Legend formatting
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.5, "lines")
  )

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/")

ggsave("umap_genomic.pdf", Umap)


##### create upset plot ##### 

# 1. Load TCGA mutation presence/absence file
# =========
tcga <- read_delim(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/upsetTCGA.txt",
  delim = "\t"
)

# =========
# 2. Convert Yes/No to binary (0/1)
# =========
tcga <- tcga %>%
  mutate(across(c(APC, KRAS, SMAD4, TP53), ~ ifelse(. == "Yes", 1, 0)))

# =========
# 3. Create mutation combination label (A, AK, AKP, AKPS, etc.)
# =========
tcga <- tcga %>%
  rowwise() %>%
  mutate(group = paste0(
    ifelse(APC == 1, "A", ""),
    ifelse(KRAS == 1, "K", ""),
    ifelse(TP53 == 1, "P", ""),
    ifelse(SMAD4 == 1, "S", "")
  )) %>%
  ungroup()

# =========
# 4. Prepare binary mutation matrix for UpSet
# =========
mutation_matrix <- tcga %>%
  dplyr::select(APC, KRAS, TP53, SMAD4) %>% as.data.frame()


mutation_df <- as.data.frame(mutation_matrix)

# Convert numeric 0/1 to logical
mutation_df[] <- lapply(mutation_df, function(x) as.logical(x))
# =========
# 5. Generate UpSet plot (UpSetR syntax)
# =========
  ComplexUpset::upset(mutation_df,
    intersect = c("APC", "KRAS", "TP53", "SMAD4"),
    base_annotations = list(
      "Intersection size" = intersection_size(
        aes(fill = after_stat())   # <- this works
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "APC"   = "steelblue",
      "KRAS"  = "tomato",
      "TP53"  = "seagreen3",
      "SMAD4" = "gold3"
    )
  )

png("upset_full.png", width = 3000, height = 2000, res = 300)
print(p)
dev.off()


###now prepare this data for ML unsupervised learning ### 
#the data is already batch corrected - so key is going to be converting mouse genes to human hugo ids so that they match 

map_mouse_to_human_expr <- function(mouse_genes, human_expr = NULL, mirror = "uswest", version = NULL) {
  # --- Connect to Ensembl marts ---g
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  # --- Step 1: get basic mouse -> human mapping (IDs only) ---
  mapping <- getLDS(
    attributes = c("mgi_symbol"),                 # mouse attributes
    filters = "mgi_symbol",
    values = mouse_genes,
    mart = mouse,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),  # human attributes
    martL = human
  )
  colnames(mapping) <- c("mouse_symbol", "human_ensembl", "human_symbol")
  
  # --- Step 2: pull human orthology type (Homology page) ---
  orthology_info <- getBM(
    attributes = c("ensembl_gene_id", "mmusculus_homolog_orthology_type"),
    mart = human
  )
  colnames(orthology_info) <- c("human_ensembl", "orthology_type")
  
  # --- Step 3: pull human Entrez IDs (Features page) ---
  id_map <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    mart = human
  )
  colnames(id_map) <- c("human_ensembl", "human_entrez")
  
  # --- Step 4: join all tables ---
  mapping_full <- mapping %>%
    left_join(orthology_info, by = "human_ensembl") %>%
    left_join(id_map, by = "human_ensembl") %>%
    filter(orthology_type == "ortholog_one2one", !is.na(human_entrez)) %>%
    distinct(mouse_symbol, .keep_all = TRUE)
  
  return(mapping_full)
  
}
mouse_genes <- (read.xlsx("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/top_genes.xlsx"))$top_genes

mapping_1to1 <- map_mouse_to_human_expr(mouse_genes)

##### now annotate human data frame 
batch_corrected_combined_mat_filt <- as.data.frame(batch_corrected_combined_mat_filt)


batch_corrected_combined_mat_filt$from_mouse <- rownames(batch_corrected_combined_mat_filt) %in% as.character(mapping_1to1$human_entrez)

batch_corrected_combined_mat_filt_small <- batch_corrected_combined_mat_filt[batch_corrected_combined_mat_filt$from_mouse == TRUE, ]


batch_corrected_combined_mat_filt_small <- batch_corrected_combined_mat_filt_small %>% rownames_to_column("human_entrez")

mapping_1to1_small <- mapping_1to1[,c("human_entrez", "human_symbol")]

batch_corrected_combined_mat_filt_small <- merge(batch_corrected_combined_mat_filt_small,mapping_1to1_small, by= "human_entrez", all.x=TRUE)

batch_corrected_combined_mat_filt_small <- batch_corrected_combined_mat_filt_small %>% column_to_rownames("human_symbol")


batch_corrected_combined_mat_filt_small <- batch_corrected_combined_mat_filt_small %>% dplyr::select(
  -from_mouse, -human_entrez
)


## now make umap plot based on this new smaller dataframe 


umap_res_small <- umap(t(batch_corrected_combined_mat_filt_small))

umap_df_small <- as.data.frame(umap_res_small$layout)

umap_df_small$sample <- colnames(batch_corrected_combined_mat_filt_small)

umap_df_small <- umap_df_small %>%
  left_join(sample_to_group, by = "sample")

Umap_small <- ggplot(umap_df_small, aes(V1, V2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  labs(color = "Mutation Group")

Umap_small <- Umap_small +
  scale_color_manual(
    values = c(
      "AK"   = "#E41A1C",
      "AKP"  = "#377EB8",
      "AKPS" = "#4DAF4A",
      "Other"= "grey70"
    )
  )

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/")

ggsave("umap_small_n=500.pdf", Umap_small)



write.xlsx(batch_corrected_combined_mat_filt, "finalMLTrainingSet.xlsx", row.names=TRUE)




####try writing new pipeline for only TCGA data ##### 
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(randomForest)
library(umap)
library(readxl)
library(openxlsx)

expr <- read.delim(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_tcga_gdc/data_mrna_seq_tpm.txt",
  row.names = 1, check.names = FALSE
)


expr <- expr %>%
  rownames_to_column("Entrez_Gene_Id") %>%
  filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Entrez_Gene_Id")

expr_log <- log2(expr + 1)

mouse_genes <- (read.xlsx("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/top_genes.xlsx"))$top_genes

mapping_1to1 <- map_mouse_to_human_expr(mouse_genes)

expr_log <- as.data.frame(expr_log)

expr_log$from_mouse <- rownames(expr_log) %in% as.character(mapping_1to1$human_entrez)


expr_log_small <- expr_log[expr_log$from_mouse == TRUE, ]

expr_log_small <- expr_log_small %>% rownames_to_column("human_entrez")

mapping_1to1_small <- mapping_1to1[,c("human_entrez", "human_symbol")]

expr_log_small <- merge(expr_log_small,mapping_1to1_small, by= "human_entrez", all.x=TRUE)

expr_log_small <- expr_log_small %>% column_to_rownames("human_symbol")

expr_log_small <- expr_log_small %>% dplyr::select(
  -from_mouse, -human_entrez
)



gene_var <- apply(expr_log_small, 1, var)
# Keep top 300 most variable genes (adjust as needed)
top_genes_var <- names(sort(gene_var, decreasing = TRUE)[1:min(300, length(gene_var))])
expr_log_small <- expr_log_small[top_genes_var, ]

expr_t <- as.data.frame(t(expr_log_small))

expr_t$sample <- substr(rownames(expr_t),1, 15)

expr_t <- expr_t %>% left_join(sample_to_group, by = "sample")

colnames(expr_t) <- make.names(colnames(expr_t))

expr_t$group <- factor(expr_t$group)
rf_fit <- randomForest(group ~ ., data = expr_t %>% dplyr::select(-sample, -NA.), 
                       ntree = 1000, importance = TRUE, proximity = TRUE)

rf_pred <- predict(rf_fit, type = "prob")  # samples x groups
rownames(rf_pred) <- expr_t$sample
posteriors <- as.data.frame(rf_pred)
#PCA on probabilities 

pca_res <- prcomp(posteriors, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[,1:2])
pca_df$sample <- rownames(posteriors)
pca_df <- pca_df %>% left_join(sample_to_group, by = "sample")

PCA_plot <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(color = "Mutation Group", title = "Supervised PCA on Random Forest Probabilities")

PCA_plot <- PCA_plot + 
  scale_color_manual(
    values = c(
      "AK"   = "#E41A1C",
      "AKP"  = "#377EB8",
      "AKPS" = "#4DAF4A",
      "Other"= "grey70"
    )
  )

# --- Step 9B: UMAP on probabilities ---
umap_res <- umap(posteriors)
umap_df <- as.data.frame(umap_res$layout)
umap_df$sample <- rownames(posteriors)
umap_df <- umap_df %>% left_join(sample_to_group, by = "sample")

UMAP_plot <- ggplot(umap_df, aes(V1, V2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(color = "Mutation Group", title = "Supervised UMAP on Random Forest Probabilities")

UMAP_plot <- UMAP_plot + scale_color_manual(
  values = c(
    "AK"   = "#E41A1C",
    "AKP"  = "#377EB8",
    "AKPS" = "#4DAF4A",
    "Other"= "grey70"
  )
)

# --- Step 10: Save plots and probabilities ---
setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/")

ggsave("RF_Supervised_PCA_n=500.pdf", PCA_plot, width = 6, height = 5)
ggsave("RF_Supervised_UMAP_n=500.pdf", UMAP_plot, width = 6, height = 5)

write.xlsx(posteriors, "RF_posterior_probabilities.xlsx", row.names = TRUE)


############# now work on correlation matrix ############ 

res_AK_small <- read.xlsx("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AK_small.xlsx")
res_AKP_small <- read.xlsx("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AKP_small.xlsx")
res_AKPS_small <- read.xlsx("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AKPS_small.xlsx")

colnames(res_AK_small) <- c("mouse_symbol", "weight")
colnames(res_AKP_small) <- c("mouse_symbol", "weight")
colnames(res_AKPS_small) <- c("mouse_symbol", "weight")

map_mouse_to_human_expr_matrix <- function(mouse_df, 
                                    keep_entrez = TRUE,
                                    aggregate = TRUE) {
  # mouse_df must contain at least:
  #   mouse_symbol : MGI gene symbol
  #   weight       : numeric score (e.g. log2FC)
  
  stopifnot(all(c("mouse_symbol", "weight") %in% colnames(mouse_df)))
  
  library(biomaRt)
  library(dplyr)
  
  # --- Connect to Ensembl marts ---
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/")
  
  # --- Step 1: mouse → human mapping ---
  mapping <- getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = unique(mouse_df$mouse_symbol),
    mart = mouse,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),
    martL = human
  )
  colnames(mapping) <- c("mouse_symbol", "human_ensembl", "human_symbol")
  
  # --- Step 2: orthology type ---
  orthology_info <- getBM(
    attributes = c("ensembl_gene_id", "mmusculus_homolog_orthology_type"),
    mart = human
  )
  colnames(orthology_info) <- c("human_ensembl", "orthology_type")
  
  # --- Step 3: human Entrez IDs ---
  if (keep_entrez) {
    id_map <- getBM(
      attributes = c("ensembl_gene_id", "entrezgene_id"),
      mart = human
    )
    colnames(id_map) <- c("human_ensembl", "human_entrez")
  } else {
    id_map <- data.frame(human_ensembl = character(),
                         human_entrez = integer())
  }
  
  # --- Step 4: join mappings ---
  mapping_full <- mapping %>%
    left_join(orthology_info, by = "human_ensembl") %>%
    left_join(id_map, by = "human_ensembl") %>%
    filter(orthology_type == "ortholog_one2one") %>%
    distinct(mouse_symbol, .keep_all = TRUE)
  
  # --- Step 5: join with weights ---
  output <- mouse_df %>%
    inner_join(mapping_full, by = "mouse_symbol") %>%
    dplyr::select(mouse_symbol, human_symbol, human_ensembl, human_entrez, weight, dplyr::everything())
  
  # --- Step 6: aggregate by human gene if requested ---
  if (aggregate) {
    output <- output %>%
      group_by(human_symbol, human_ensembl, human_entrez) %>%
      summarize(weight = mean(weight, na.rm = TRUE), .groups = "drop")
  }
  
  return(output)
}

 res_AK_small_complete <- map_mouse_to_human_expr_matrix(res_AK_small)
 res_AKP_small_complete <- map_mouse_to_human_expr_matrix(res_AKP_small)
 res_AKPS_small_complete <- map_mouse_to_human_expr_matrix(res_AKPS_small)
 
 ##### now use function to map each huuman sample to a group ##### 
 
 res_AK_small_complete <- res_AK_small_complete %>%
   dplyr::filter(!is.na(human_entrez)) %>%
   distinct(human_entrez, .keep_all = TRUE) %>%
   column_to_rownames("human_entrez") %>%
   dplyr::select(weight)
 
 
 res_AKP_small_complete <- res_AKP_small_complete %>%
   dplyr::filter(!is.na(human_entrez)) %>%
   distinct(human_entrez, .keep_all = TRUE) %>%
   column_to_rownames("human_entrez") %>%
   dplyr::select(weight)
 
 res_AKPS_small_complete <- res_AKPS_small_complete %>%
   dplyr::filter(!is.na(human_entrez)) %>%
   distinct(human_entrez, .keep_all = TRUE) %>%
   column_to_rownames("human_entrez") %>%
   dplyr::select(weight)
 
 AK_weights <-  res_AK_small_complete$weight
 names(AK_weights) <- rownames(res_AK_small_complete)
 
 AKP_weights <-  res_AKP_small_complete$weight
 names(AKP_weights) <- rownames(res_AKP_small_complete)
 
 AKPS_weights <-  res_AKPS_small_complete$weight
 names(AKPS_weights) <- rownames(res_AKPS_small_complete)
 
 cor_score_per_sample <- function(weights, expr_mat) {
   # take intersection
   common <- intersect(names(weights), rownames(expr_mat))
   if(length(common) < 5) stop("Too few genes overlapping.")
   w <- weights[common]
   e <- expr_mat[common, , drop = FALSE]   # rows genes, cols samples
   # compute Pearson correlation between weight vector and each sample's expression vector
   apply(e, 2, function(sample_vec) cor(as.numeric(w), as.numeric(sample_vec), use = "pairwise.complete.obs", method="spearman"))
 }
 
 scoreAK_cor_TCGA   <- cor_score_per_sample(AK_weights, TCGA_zscores)
 scoreAKP_cor_TCGA   <- cor_score_per_sample(AKP_weights, TCGA_zscores)
 scoreAKPS_cor_TCGA   <- cor_score_per_sample(AKPS_weights, TCGA_zscores)
 
 scores_cor_df <- data.frame(
   sample_id = colnames(TCGA_zscores),
   AK   = scoreAK_cor_TCGA,
   AKP  = scoreAKP_cor_TCGA,
   AKPS = scoreAKPS_cor_TCGA,
   stringsAsFactors = FALSE
 )
 scores_cor_df$pred_label <- apply(scores_cor_df[, c("AK","AKP","AKPS")], 1, function(x) c("AK","AKP","AKPS")[which.max(x)])
 
 scores_cor_df$maxCorr <- apply(scores_cor_df[,c("AK","AKP","AKPS")], 1, max)
 
 
 write.xlsx(scores_cor_df, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/corrMatrix_n=1000.xlsx")
 
 
 ###########
 
 
 