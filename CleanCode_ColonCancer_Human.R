
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(tidyr)
library(ggplot2)
source("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Code/Helperfunctions_Human.R")
library(CMSclassifier)
library(DESeq2)
library(readr)
library(ComplexUpset)
library(ComplexHeatmap)
library(circlize)
library(umap)
library(grid)


# ---------------------------------------------------------------
# 1. Download raw counts
# ---------------------------------------------------------------
message("1. Downloading TCGA-COAD data...")

coad_query <- GDCquery(
  project       = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query = coad_query)
coad_data <- GDCprepare(query = coad_query, summarizedExperiment = TRUE)

# ---------------------------------------------------------------
# 2. Gene annotations
# ---------------------------------------------------------------
message("2. Extracting gene annotations...")

gene_meta <- as.data.frame(rowRanges(coad_data)) %>%
  dplyr::select(gene_id, gene_name, gene_type, width) %>%
  mutate(width_kb = width / 1000)

# ---------------------------------------------------------------
# 3. Convert raw counts to log2(TPM + 1) -- +1 to guard against log(0)

# ---------------------------------------------------------------
message("3. Converting raw counts to log2(TPM + 1)...")

counts_long <- as.data.frame(assay(coad_data)) %>%
  rownames_to_column("ensembl_id") %>%
  left_join(
    gene_meta %>% dplyr::select(gene_id, width_kb),
    by = c("ensembl_id" = "gene_id")
  ) %>%
  dplyr::filter(!is.na(width_kb), width_kb > 0)  # need a valid length to compute TPM

sample_cols <- setdiff(names(counts_long), c("ensembl_id", "width_kb"))

coad_tpm <- counts_to_tpm(counts_long, sample_cols) %>%
  dplyr::select(ensembl_id, all_of(sample_cols)) %>%
  column_to_rownames("ensembl_id") %>%
  as.matrix()

coad_tpm_log2 <- log2(coad_tpm + 1)

# ---------------------------------------------------------------
# 4. Keep one sample per patient
#    When a patient has multiple aliquots, keep the one with the
#    lexicographically-highest portion/analyte code.
# ---------------------------------------------------------------
message("4. Selecting one sample per patient...")

sample_info <- tibble(
  full_id = colnames(coad_tpm_log2),
  base_id = substr(full_id, 1, 15),
  # portion/analyte block: the 5th hyphen-delimited barcode field, e.g. "11R"
  portion = sapply(strsplit(full_id, "-"), function(x) if (length(x) >= 5) x[5] else NA_character_)
)

n_before <- nrow(sample_info)

sample_info_unique <- sample_info %>%
  dplyr::filter(!is.na(portion)) %>%
  group_by(base_id) %>%
  slice_max(order_by = portion, n = 1, with_ties = FALSE) %>%
  ungroup()

message(sprintf(
  "   %d aliquots -> %d unique patients (%d dropped: unparsed barcode)",
  n_before, nrow(sample_info_unique), n_before - nrow(sample_info_unique)
))

coad_tpm_log2_unique <- coad_tpm_log2[, sample_info_unique$full_id]
colnames(coad_tpm_log2_unique) <- sample_info_unique$base_id

write.csv(coad_tpm_log2_unique, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Main Figures/Figure 1/coad_tpm_log2_unique.csv")

# ---------------------------------------------------------------
# 5. Map Ensembl IDs to HGNC symbols for CIBERSORT
# ---------------------------------------------------------------
message("5. Mapping Ensembl IDs to HGNC symbols...")

coad_tpm_cibersort <- map_and_collapse_ids(
  coad_tpm_log2_unique,
  to_keytype = "SYMBOL",
  label      = "HGNC"
)

message("Done. CIBERSORT input prepared with HGNC symbols as row names: ",
        nrow(coad_tpm_cibersort), " genes x ", ncol(coad_tpm_cibersort), " samples.")


###### now perform statistical analysis CIBERSORT 

# Load sample list
A_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/A.csv")
AK_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AK_genomic.csv")
AKP_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AKP_genomic.csv")
AKPS_TCGA <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Genomic/AKPS_genomic.csv")



A_samples <- A_TCGA$Sample.ID
AK_samples <- AK_TCGA$Sample.ID
AKP_samples <- AKP_TCGA$Sample.ID
AKPS_samples <- AKPS_TCGA$Sample.ID


# Now, subset this HGNC-symbol-mapped TPM data into AK and AKP groups
COAD_counts_tpm_AK_for_cibersort <- coad_tpm_cibersort[,colnames(coad_tpm_cibersort) %in% AK_samples]
COAD_counts_tpm_AKP_for_cibersort <- coad_tpm_cibersort[,colnames(coad_tpm_cibersort) %in% AKP_samples]
COAD_counts_tpm_A_for_cibersort <- coad_tpm_cibersort[,colnames(coad_tpm_cibersort) %in% A_samples]
COAD_counts_tpm_AKPS_for_cibersort <- coad_tpm_cibersort[,colnames(coad_tpm_cibersort) %in% AKPS_samples]

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
                     return(list(p.value = NA))
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
    mean_AKP = mean1,    
    mean_AKPS = mean2,   
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

########### make plot 1 A,B,C,D ############-- all using this log(TPM+1) dataframe 
#Make barplot of immunecells 
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

# Stacked bar plot
humanStackedImmune <- ggplot(df_summary,
                       aes(x = Group,
                           y = Fraction,
                           fill = CellType)) +
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7) +
  labs(
    title = "Immune Cell Fractions Across Tumor Progression",
    x = "Progression Stage",
    y = "Mean Immune Cell Fraction",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Human_immune_barplot.pdf",humanStackedImmune , width = 14, height = 8)

############# NOW PERFORM CMS SUBTYPING so we can make 1B ############### 

# ---------------------------------------------------------------
# 1. Build CMS input: log2(TPM+1) with Entrez IDs as row names
# ---------------------------------------------------------------
message("Preparing CMS classification input...")

coad_cms_input <- map_and_collapse_ids(
  coad_tpm_log2_unique,
  to_keytype = "ENTREZID",
  label      = "Entrez"
)

message("CMS input prepared: ", nrow(coad_cms_input), " genes x ",
        ncol(coad_cms_input), " samples.")

# ---------------------------------------------------------------
# 2. Run CMS classification
# ---------------------------------------------------------------
message("Running CMS classification (random forest)...")

cms_results <- classifyCMS(as.matrix(coad_cms_input), method = "RF")

message("CMS classification results:")
print(cms_results$predictedCMS)


write.csv(cms_results$predictedCMS,
          "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classification_results.csv",
          row.names = TRUE)
message("Pipeline execution complete.")

############################## now that CMS classifications are done make stacked barplot ################

CMS_classifications <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classification_results.csv")

CMS_classifications <- CMS_classifications %>% rename(
  "X" = "SampleName"
)

CMS_classifications_barplot <- CMS_classifications %>% mutate(group=case_when(
  SampleName %in% A_samples ~ "A", 
  SampleName %in% AK_samples ~ "AK", 
  SampleName %in% AKP_samples ~ "AKP", 
  SampleName %in% AKPS_samples ~ "AKPS"
))

CMS_classifications_barplot <- CMS_classifications_barplot[!is.na(CMS_classifications_barplot$group),]

CMS_classifications_barplot$group <- factor(CMS_classifications_barplot$group, levels = c("A","AK","AKP","AKPS"))

write.csv(CMS_classifications_barplot, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classifications_alluvial.csv")

CMS_classifications_barplot <- CMS_classifications_barplot %>% 
  mutate(RF = replace_na(RF, "Other"))

CMS_classifications_barplot <- CMS_classifications_barplot[!is.na(CMS_classifications_barplot$group),]

CMS_classifications_barplot_grouped <- CMS_classifications_barplot %>%
  group_by(RF, group) %>%
  summarise(Freq =n(), .groups = "drop")

CMS_classifications_barplot_grouped <- CMS_classifications_barplot_grouped %>% group_by(group) %>% mutate(
  proportion = Freq/sum(Freq) 
) %>% ungroup()


# Stacked bar plot
humanStackedCMS <- ggplot(
  CMS_classifications_barplot_grouped,
  aes(x = group, y = proportion, fill = RF)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  labs(
    title = "CMS Subtype Composition Across Tumor Progression",
    x = "Progression Stage",
    y = "Fraction of Samples Belonging to CMS Subtype",
    fill = "CMS Subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

humanStackedCMS


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Human_barplot_CMS.pdf",humanStackedCMS, width = 14, height = 8)

################ now make dds object ################
#===============================================================
# TCGA-COAD: Build a DESeqDataSet with A / AK / AKP / AKPS groups
#===============================================================

# ---------------------------------------------------------------
# 1. Assemble the grouping variable from the four ID lists
# ---------------------------------------------------------------
message("Building A/AK/AKP/AKPS grouping variable...")

group_lookup <- bind_rows(
  tibble(base_id = A_samples,    group = "A"),
  tibble(base_id = AK_samples,   group = "AK"),
  tibble(base_id = AKP_samples,  group = "AKP"),
  tibble(base_id = AKPS_samples, group = "AKPS")
)

# Sanity checks: no sample listed in more than one group, and every
# ID is exactly 15 characters as expected
dup_ids <- group_lookup$base_id[duplicated(group_lookup$base_id)]
if (length(dup_ids) > 0) {
  stop("These sample IDs appear in more than one group list: ",
       paste(unique(dup_ids), collapse = ", "))
}

bad_length <- group_lookup$base_id[nchar(group_lookup$base_id) != 15]
if (length(bad_length) > 0) {
  warning("These IDs are not 15 characters long: ",
          paste(bad_length, collapse = ", "))
}

# ---------------------------------------------------------------
# 2. Keep every de-duplicated sample; anything not in one of the
#    four lists is assigned to "Other"
# ---------------------------------------------------------------
message("Building raw count matrix for DESeq2...")

sample_info_grouped <- sample_info_unique %>%
  left_join(group_lookup, by = "base_id") %>%
  mutate(group = if_else(is.na(group), "Other", group))

n_other <- sum(sample_info_grouped$group == "Other")
if (n_other > 0) {
  message(sprintf(
    "%d de-duplicated samples are not in any of the A/AK/AKP/AKPS lists and were assigned to \"Other\".",
    n_other
  ))
}

count_mat <- counts_long %>%
  dplyr::select(ensembl_id, all_of(sample_info_grouped$full_id)) %>%
  column_to_rownames("ensembl_id") %>%
  as.matrix()

storage.mode(count_mat) <- "integer"  # DESeq2 requires integer counts

colnames(count_mat) <- sample_info_grouped$base_id[
  match(colnames(count_mat), sample_info_grouped$full_id)
]

# ---------------------------------------------------------------
# 3. Sample metadata (colData)
# ---------------------------------------------------------------
col_data <- sample_info_grouped %>%
  dplyr::select(base_id, group) %>%
  mutate(group = factor(group, levels = c("A", "AK", "AKP", "AKPS", "Other"))) %>%
  column_to_rownames("base_id")

col_data <- col_data[colnames(count_mat), , drop = FALSE]  # match order exactly

stopifnot(identical(rownames(col_data), colnames(count_mat)))

# ---------------------------------------------------------------
# 4. Build the DESeqDataSet
# ---------------------------------------------------------------
message("Creating DESeqDataSet...")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = col_data,
  design    = ~ group
)

message("DESeqDataSet created: ", nrow(dds), " genes x ", ncol(dds), " samples.")
print(table(col_data$group))

vsd <- vst(dds)

data_human <- as.data.frame(assay(vsd))

######## now add clinical data to address reviewer 3 #########

tcga_clinical <- as.data.frame(colData(coad_data))
tcga_clinical$portion <- sapply(strsplit(rownames(tcga_clinical), "-"), function(x) if (length(x) >= 5) x[5] else NA_character_)
tcga_clinical$sample <- substr(tcga_clinical$sample,1,15)

tcga_clinical <- tcga_clinical%>%
  group_by(sample) %>%
  slice_max(order_by = portion, n = 1, with_ties = FALSE) %>%
  ungroup()

clinical <- tcga_clinical %>%
  dplyr::select(
    sample,
    ajcc_pathologic_stage,
    paper_MSI_status
  )


# Run UMAP on samples (transpose so rows = samples)
set.seed(123)
umap_result_human <- uwot::umap(
  t(data_human),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "euclidean"
)

# Convert to data frame
umap_df_human <- as.data.frame(umap_result_human)
colnames(umap_df_human) <- c("UMAP1", "UMAP2")

umap_df_human <- umap_df_human %>% rownames_to_column("sample")

umap_df_human <- umap_df_human <- umap_df_human %>%
  left_join(clinical, by="sample")


# Assign groups
umap_df_human <- umap_df_human %>%
  mutate(grp_gene = case_when(
    sample %in% A_samples ~ "A",
    sample %in% AK_samples ~ "AK",
    sample %in% AKP_samples ~ "AKP",
    sample %in% AKPS_samples ~ "AKPS"
  ))


umap_df_human <- umap_df_human %>%
  mutate(
    # Collapse AJCC stage
    stage_group = case_when(
      grepl("^Stage IV", ajcc_pathologic_stage) ~ "Stage IV",
      grepl("^Stage III", ajcc_pathologic_stage) ~ "Stage III",
      grepl("^Stage II", ajcc_pathologic_stage) ~ "Stage II",
      grepl("^Stage I", ajcc_pathologic_stage) ~ "Stage I",
      TRUE ~ NA_character_
    ),
    
    # Clean MSI status
    MSI_group = case_when(
      paper_MSI_status %in% c("MSS", "MSI-H", "MSI-L") ~ paper_MSI_status,
      TRUE ~ NA_character_
    )
  )

umap_df_human <- umap_df_human %>%
  mutate(
    MSI_binary = case_when(
      paper_MSI_status == "MSI-H" ~ "MSI-H",
      paper_MSI_status %in% c("MSS", "MSI-L") ~ "MSS/MSI-L",
      TRUE ~ NA_character_
    )
  )


# Plot
p_genotype <- ggplot(umap_df_human,
                     aes(x = UMAP1, y = UMAP2, color = grp_gene)) +
  geom_point(size = 7) +
  labs(
    title = "UMAP Plot",
    x = "UMAP1",
    y = "UMAP2",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.ticks = element_line(linewidth = 0.5)
  )


p_msi <- ggplot(umap_df_human,
                     aes(x = UMAP1, y = UMAP2, color = MSI_binary)) +
  geom_point(size = 7) +
  labs(
    title = "UMAP Plot",
    x = "UMAP1",
    y = "UMAP2",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.ticks = element_line(linewidth = 0.5)
  )

p_stage <- ggplot(umap_df_human,
                aes(x = UMAP1, y = UMAP2, color = stage_group)) +
  geom_point(size = 7) +
  labs(
    title = "UMAP Plot",
    x = "UMAP1",
    y = "UMAP2",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.ticks = element_line(linewidth = 0.5)
  )

library(patchwork)

umapCombined <- (p_genotype | p_msi | p_stage)

ggsave(
  filename = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/updatedUmap.pdf", 
  plot = umapCombined,                 
  width = 15, 
  height = 6, 
  units = "in", 
  dpi = 300
)



write.csv(umap_df_human, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 1/Figure B/umap.csv")

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/umap.pdf", umap_human, width = 14, height = 8)

######## now make upset plot ############ 
tcga <- read_delim(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 1/Figure A/upsetTCGA.txt",
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
p <- ComplexUpset::upset(
  mutation_df,
  intersect = c("APC","KRAS","TP53","SMAD4"),
  base_annotations = list(
    "Intersection size" = intersection_size()
  )
)

png("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/upsetPlot.png",
    width = 3000, height = 2000, res = 300)

print(p)

dev.off()

############# now make figure 1D which will be a heatmap showing pairwise correlations ######### -- Using VST transformed data (data_human)

rv <- matrixStats::rowVars(as.matrix(data_human))

data_human_variable <- data_human[order(rv, decreasing = TRUE)[1:5000], ]

cor_mat <- cor(data_human_variable, method="spearman") 
#meta data is col_data 
stopifnot(identical(rownames(col_data), colnames(data_human_variable)))

write.csv(data_human, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 1/HumanData.csv")
write.csv(cor_mat, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 1/Figure E/corMatrix.csv")

group_colors <- c(
  A     = "#377EB8",   # blue
  AK    = "#4DAF4A",   # green
  AKP   = "#000000",   # black
  AKPS  = "#E41A1C",   # red
  Other = "grey80"
)

ha <- HeatmapAnnotation(
  Genotype = col_data$group,
  col = list(Genotype = group_colors),
  annotation_name_side = "left"
)

hc <- hclust(as.dist(1 - cor_mat), method = "complete")
###########################################################
## Heatmap
###########################################################

ht <- Heatmap(
  cor_mat,
  
  name = "Spearman\nρ",
  
  top_annotation = ha,
  
  cluster_rows = hc,
  cluster_columns = hc,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  
  col = colorRamp2(
    c(0.5, 0.75, 1),
    c("white", "orange", "darkred")
  ),
  
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation"
  )
)

png("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/humanHeatMap.png",
    width = 3000, height = 2000, res = 300)
draw(ht)

dev.off()


############################ now see if predicted transcriptomic classes (AK-like, AKP-like, AKPS-like) separate the data better than genomic labels #########

output_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision"
fig8_dir   <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Datasets"
gsea_human_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Transcriptomic Analysis"
tcga_zscore_file <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"

for (d in c(output_dir, fig8_dir, gsea_human_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)


pred_results <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/pred_results.csv")

pred_results_15 <- pred_results %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE)  # guards against duplicate-key row multiplication in any join below

dup_after_trunc <- sum(duplicated(substr(pred_results$Sample, 1, 15)))
if (dup_after_trunc > 0) {
  warning(sprintf("%d pred_results samples collapsed to a duplicate 15-char ID after truncation; kept first occurrence only.",
                  dup_after_trunc))
}

########################################################
# 1. UMAP: predicted transcriptomic class vs. sample separation
########################################################
message("1. Building UMAP...")


tcga_zscores <- read_delim(tcga_zscore_file, delim = "\t", show_col_types = FALSE) %>%
  filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

TCGA_zscores_t <- t(tcga_zscores) %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

TCGA_zscores_t$Sample <- substr(TCGA_zscores_t$Sample, 1, 15)

common_samples <- intersect(TCGA_zscores_t$Sample, pred_results_15$Sample)
TCGA_zscores_t <- TCGA_zscores_t[TCGA_zscores_t$Sample %in% common_samples, ]

pca_input <- left_join(TCGA_zscores_t, pred_results_15[, c("Sample", "Predicted_Class")], by = "Sample")

pca_input <- pca_input %>% column_to_rownames("Sample")

numeric_data <- pca_input %>% dplyr::select(-Predicted_Class)

numeric_data <- numeric_data[
  , colSums(is.na(numeric_data)) < nrow(numeric_data),
  drop = FALSE
]

set.seed(123)
umap_res <- umap(numeric_data)



umap_df <- data.frame(Sample = rownames(numeric_data), UMAP1 = umap_res$layout[, 1], UMAP2 = umap_res$layout[, 2]) %>%
  left_join(pred_results_15[, c("Sample", "Predicted_Class")], by = "Sample") %>%
  mutate(Predicted_Class = factor(Predicted_Class, levels = c("AK", "AKP", "AKPS", "None")))

class_colors <- c(AK = "#E41A1C", AKP = "#377EB8", AKPS = "#4DAF4A", "None" = "grey70")

UMAP_Transcriptomic <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Predicted_Class)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Predicted_Class), type = "norm", level = 0.95, geom = "polygon", alpha = 0.15, color = NA) +
  scale_color_manual(values = class_colors) +
  scale_fill_manual(values = class_colors) +
  theme_minimal() +
  labs(title = "", x = "UMAP1", y = "UMAP2", color = "Mutation Group", fill = "Mutation Group") +
  theme(legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14), legend.key.size = unit(1.5, "lines"))

ggsave(file.path(output_dir, "umap_transcriptomic_Aug27.pdf"), UMAP_Transcriptomic, width = 8, height = 8)

write.csv(TCGA_zscores_t, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/TCGA_zscores_t.csv")
write.csv(umap_df , "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/umap_df.csv")
########################################################
# 3. QC heatmap: module eigengenes ordered by predicted class
########################################################
message("3. Building module-eigengene QC heatmap...")

eigengenes_annot <- eigengenes_human_mat %>%
  as.data.frame() %>%
  mutate(Sample = substr(rownames(.), 1, 15)) %>%   # normalize to 15-char base ID, matching pred_results_15
  left_join(pred_results_15, by = "Sample") %>%
  mutate(Predicted_Class = factor(Predicted_Class, levels = c("AK", "AKP", "AKPS", "None"))) %>%
  arrange(Predicted_Class) %>%
  dplyr::filter(Predicted_Class != "None")

mat_ordered <- eigengenes_annot %>% dplyr::select(starts_with("Module")) %>% as.matrix()
rownames(mat_ordered) <- eigengenes_annot$Sample

row_anno <- data.frame(Predicted_Class = eigengenes_annot$Predicted_Class, row.names = eigengenes_annot$Sample)
ann_colors <- list(Predicted_Class = c(AK = "#1b9e77", AKP = "#d95f02", AKPS = "#7570b3"))

pheatmap(
  mat_ordered, cluster_rows = FALSE, cluster_cols = TRUE,
  annotation_row = row_anno, annotation_colors = ann_colors,
  scale = "row", show_rownames = FALSE, border_color = NA,
  fontsize_col = 11, main = "Module Expression by Predicted Class",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  filename = file.path(output_dir, "module_eigengene_QC_heatmap.png"),
  width = 8, height = 6
)

write.csv(mat_ordered, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/Figure 7B/mat_ordered.csv")

########################################################
# 4. Confusion matrix: genomic vs. transcriptomic classification
########################################################
message("4. Building confusion matrix...")

genomic_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data"
transcriptomic_dir <- file.path(genomic_dir, "CIBERSORT/Oct 12/lists/Transcriptomic")

AK_samples_genomic   <- read.csv(file.path(genomic_dir, "AK.csv"))$Sample.ID
AKP_samples_genomic  <- read.csv(file.path(genomic_dir, "AKP.csv"))$Sample.ID
AKPS_samples_genomic <- read.csv(file.path(genomic_dir, "AKPS.csv"))$Sample.ID



AK_samples_transcriptomic   <- substr(pred_results[pred_results$Predicted_Class=="AK","Sample"],1,15)
AKP_samples_transcriptomic  <- substr(pred_results[pred_results$Predicted_Class=="AKP","Sample"],1,15)
AKPS_samples_transcriptomic <- substr(pred_results[pred_results$Predicted_Class=="AKPS","Sample"],1,15)
None_samples_transcriptomic <- substr(pred_results[pred_results$Predicted_Class=="None","Sample"],1,15)

genomic_df <- data.frame(
  Sample = c(AK_samples_genomic, AKP_samples_genomic, AKPS_samples_genomic),
  Genomic_Class = c(rep("AK", length(AK_samples_genomic)), rep("AKP", length(AKP_samples_genomic)),
                    rep("AKPS", length(AKPS_samples_genomic)))
)
transcriptomic_df <- data.frame(
  Sample = c(AK_samples_transcriptomic, AKP_samples_transcriptomic, AKPS_samples_transcriptomic, None_samples_transcriptomic),
  Transcriptomic_Class = c(rep("AK", length(AK_samples_transcriptomic)), rep("AKP", length(AKP_samples_transcriptomic)),
                           rep("AKPS", length(AKPS_samples_transcriptomic)), rep("None", length(None_samples_transcriptomic)))
)

class_compare <- merge(genomic_df, transcriptomic_df, by = "Sample", all = TRUE) %>%
  dplyr::filter(Sample %in% genomic_df$Sample)

table_compare <- table(class_compare$Genomic_Class, class_compare$Transcriptomic_Class)

expected_cols <- c("None", "AK", "AKP", "AKPS")
missing_cols <- setdiff(expected_cols, colnames(table_compare))
for (mc in missing_cols) table_compare <- cbind(table_compare, setNames(data.frame(rep(0, nrow(table_compare))), mc))
table_compare <- table_compare[, expected_cols]


png(file.path(output_dir, "ConfusionMatrix_Aug28.png"), width = 6, height = 6, units = "in", res = 300)
ph <- pheatmap(
  table_compare, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
  color = colorRampPalette(c("white", "steelblue"))(50), fontsize = 17, fontsize_number = 15,
  angle_col = 0, number_format = "%.0f", cellwidth = 60, cellheight = 60, silent = TRUE
)
grid.newpage()
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.draw(ph$gtable)
grid.text("Transcriptomic Classification", x = 0.5, y = unit(0.02, "npc"), gp = gpar(fontsize = 11))
grid.text("Genomic Classification", x = unit(0.02, "npc"), y = 0.5, rot = 90, gp = gpar(fontsize = 11))
dev.off()

write.csv(table_compare, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/Figure 7C/tableCompare.csv")


########################################################
# 5. GSEA on transcriptomic (predicted) labels: AKP vs AK, AKPS vs AKP
########################################################
message("5. Preparing GSEA input from transcriptomic labels...")

raw_counts <- assay(coad_data)
.
col_info <- tibble(
  full_id = colnames(raw_counts),
  base_id = substr(full_id, 1, 15),
  portion = sapply(strsplit(full_id, "-"), function(x) if (length(x) >= 5) x[5] else NA_character_)
) %>%
  dplyr::filter(!is.na(portion)) %>%
  group_by(base_id) %>%
  slice_max(order_by = portion, n = 1, with_ties = FALSE) %>%
  ungroup()

raw_counts_unique <- raw_counts[, col_info$full_id]
colnames(raw_counts_unique) <- col_info$base_id

human_sample_metadata <- data.frame(
  sample_id = colnames(raw_counts_unique),
  group = case_when(
    colnames(raw_counts_unique) %in% AK_samples_transcriptomic  ~ "AK",
    colnames(raw_counts_unique) %in% AKP_samples_transcriptomic ~ "AKP",
    colnames(raw_counts_unique) %in% AKPS_samples_transcriptomic ~ "AKPS",
    TRUE ~ NA_character_
  )
) %>% dplyr::filter(!is.na(group))

dds_human_transcriptomic <- DESeqDataSetFromMatrix(
  countData = raw_counts_unique[, human_sample_metadata$sample_id],
  colData   = human_sample_metadata,
  design    = ~ group
)
stopifnot(identical(colnames(dds_human_transcriptomic), human_sample_metadata$sample_id))  # cheap safety net
dds_human_transcriptomic <- DESeq(dds_human_transcriptomic)

human_gsea_comparisons <- list(c("AKP", "AK"), c("AKPS", "AKP"))

for (cmp in human_gsea_comparisons) {
  res_df <- results(dds_human_transcriptomic, contrast = c("group", cmp[1], cmp[2])) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_meta %>% dplyr::select(gene_id, gene_name), by = c("ensembl_id" = "gene_id")) %>%
    dplyr::filter(!is.na(stat) & !is.na(gene_name)) %>%
    arrange(desc(stat)) %>%    # NOTE: ranked by Wald stat here, vs. sign(LFC)*-log10(p) in the mouse
    dplyr::select(gene_name, stat)  # pipeline -- both are legitimate but not directly comparable in magnitude
  
  fname <- file.path(gsea_human_dir, paste0(cmp[1], "_vs_", cmp[2], "_ranked_genes2.rnk"))
  write.table(res_df, file = fname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ---------------------------------------------------------------
read_human_gsea_result <- function(label, base_dir, baseline) {
  run_dir  <- find_gsea_output_dir(base_dir, label, baseline) 
  neg_file <- list.files(run_dir, pattern = "^gsea_report_for_na_neg_.*\\.tsv$", full.names = TRUE)
  pos_file <- list.files(run_dir, pattern = "^gsea_report_for_na_pos_.*\\.tsv$", full.names = TRUE)
  bind_rows(read_tsv(neg_file, show_col_types = FALSE), read_tsv(pos_file, show_col_types = FALSE)) %>%
    dplyr::rename(pathway = NAME, padj = `FDR q-val`) %>%
    mutate(NES = as.numeric(NES), padj = as.numeric(padj)) %>%
    dplyr::filter(!is.na(NES))
}

AKvAKP    <- read_human_gsea_result("AKP", gsea_human_dir, "AK2")
AKPvAKPS  <- read_human_gsea_result("AKPS", gsea_human_dir, "AKP2")

########################################################
# 6. Hallmark pathway heatmaps (immune + oncogenic)
########################################################
message("6. Building pathway heatmaps...")

immune_pathways <- c(
  "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT", "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_TGF_BETA_SIGNALING"
)
oncogenic_pathways <- c(
  "HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_P53_PATHWAY",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

comparison_labels <- c("AKP_vs_AK" = "AKP vs. AK", "AKPS_vs_AKP" = "AKPS vs. AKP")

plot_full_hallmark_heatmap_with_priority <- function(gsea_results_list, comparisons, priority_pathways,
                                                     xlsx_path, save_path,
                                                     save_width = 18, save_height = 8,
                                                     strip_text_size = 11,
                                                     pval_cutoff = 0.05) {
  all_results <- purrr::imap_dfr(gsea_results_list, ~ {
    .x %>% dplyr::select(pathway, NES, padj) %>% mutate(comparison = .y)
  }) %>% dplyr::filter(grepl("^HALLMARK_", pathway))
  
  all_pathways <- unique(all_results$pathway)
  priority_order <- priority_pathways[priority_pathways %in% all_pathways]
  remaining <- setdiff(sort(all_pathways), priority_order)
  final_order <- c(priority_order, remaining)
  

  nes_matrix <- all_results %>%
    dplyr::select(pathway, comparison, NES) %>%
    pivot_wider(names_from = comparison, values_from = NES) %>%
    column_to_rownames("pathway")
  nes_matrix <- nes_matrix[final_order, , drop = FALSE]
  colnames(nes_matrix) <- comparisons[colnames(nes_matrix)]
  write.xlsx(nes_matrix %>% rownames_to_column("pathway"), xlsx_path)
  
  plot_data <- all_results %>%
    mutate(
      pathway = factor(pathway, levels = rev(final_order)),
      Comparison = factor(comparisons[comparison], levels = comparisons[names(gsea_results_list)]),
      Significance = if_else(padj < pval_cutoff, "*", "")  # "*" matches the mouse heatmap's symbol
    )
  
  p <- ggplot(plot_data, aes(Comparison, pathway, fill = NES)) +
    geom_tile() +
    geom_text(aes(label = Significance), color = "black", size = 9) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
    facet_wrap(~Comparison, nrow = 1, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 8), legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 8),
      strip.text = element_text(size = strip_text_size, face = "bold"),
      strip.text.x = element_text(size = strip_text_size)
    )
  
  ggsave(save_path, p, width = save_width, height = save_height)
  p
}

AKvAKP_immune      <- AKvAKP[AKvAKP$pathway %in% immune_pathways, ]
AKvAKP_oncogenic   <- AKvAKP[AKvAKP$pathway %in% oncogenic_pathways, ]
AKPvAKPS_immune    <- AKPvAKPS[AKPvAKPS$pathway %in% immune_pathways, ]
AKPvAKPS_oncogenic <- AKPvAKPS[AKPvAKPS$pathway %in% oncogenic_pathways, ]

gsea_results_list_immune <- list(AKP_vs_AK = AKvAKP_immune, AKPS_vs_AKP = AKPvAKPS_immune)
gsea_results_list_oncogenic <- list(AKP_vs_AK = AKvAKP_oncogenic, AKPS_vs_AKP = AKPvAKPS_oncogenic)

hallMarkHeatMapHuman_IMMUNE <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list_immune, comparison_labels, immune_pathways,
  xlsx_path = file.path(fig8_dir, "immuneNES.xlsx"),
  save_path = file.path(output_dir, "hallMarkHeatMapHuman_immune_TRANSCRIPTOMIC_Aug29.pdf"),
  save_width = 18, save_height = 8, pval_cutoff = 0.05
)
hallMarkHeatMapHuman_ONCOGENIC <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list_oncogenic, comparison_labels, oncogenic_pathways,  # FIX: was `immune_pathways` in the original
  xlsx_path = file.path(fig8_dir, "oncogenicNES.xlsx"),
  save_path = file.path(output_dir, "hallMarkHeatMapHuman_oncogenic_TRANSCRIPTOMIC_Aug29.pdf"),
  save_width = 18, save_height = 8, pval_cutoff = 0.05
)


########################################################
# 7. Volcano plots (transcriptomic-label groups)
########################################################
message("7. Building volcano plots...")

make_volcano_human <- function(dds, contrast, plot_title, extra_labels = NULL,
                               fc_cut = 2, p_cut = 0.05, n_label = 10, show_legend = TRUE) {
  res <- lfcShrink(dds, contrast = contrast, type = "ashr") %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_meta %>% dplyr::select(gene_id, gene_name, gene_type), by = c("ensembl_id" = "gene_id")) %>%
    dplyr::filter(gene_type == "protein_coding")
  
  rownames(res) <- make.unique(res$gene_name)
  res$fdr_adj_pvalue <- p.adjust(res$pvalue, method = "fdr")
  
  get_top_genes <- function(df, fc_filter, n = 10) {
    subset_df <- df[fc_filter(df$log2FoldChange), ]
    subset_df <- subset_df[order(subset_df$fdr_adj_pvalue), ]
    subset_df$gene_name[1:min(n, nrow(subset_df))]
  }
  up_genes   <- get_top_genes(res, function(fc) fc >= fc_cut, n_label)
  down_genes <- get_top_genes(res, function(fc) fc <= -fc_cut, n_label)
  select_labels <- unique(c(up_genes, down_genes, extra_labels, "REG4"))
  
  p <- EnhancedVolcano(
    res, lab = res$gene_name, selectLab = select_labels,
    x = "log2FoldChange", y = "fdr_adj_pvalue", title = plot_title,
    pCutoff = p_cut, FCcutoff = fc_cut, drawConnectors = TRUE,
    maxoverlapsConnectors = Inf, lengthConnectors = unit(0.0075, "npc"), boxedLabels = FALSE
  )
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}


volcano_AKPvAK   <- make_volcano_human(dds_human_transcriptomic, c("group", "AKP", "AK"),  "AKP vs AK", show_legend = TRUE)
volcano_AKPSvAKP <- make_volcano_human(dds_human_transcriptomic, c("group", "AKPS", "AKP"), "AKPS vs AKP", show_legend = FALSE)

ggsave(file.path(output_dir, "volcanoAKPvAK_TRANSCRIPTOMIC_Aug29.pdf"), volcano_AKPvAK, width = 10, height = 8)
ggsave(file.path(output_dir, "volcanoAKPSvAKP_TRANSCRIPTOMIC_Aug29.pdf"), volcano_AKPSvAKP, width = 10, height = 8)

message("Done.")

AKPvAK <- lfcShrink(dds_human_transcriptomic, contrast = c("group", "AKP", "AK"), type = "ashr") %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_meta %>% dplyr::select(gene_id, gene_name, gene_type), by = c("ensembl_id" = "gene_id")) %>%
  dplyr::filter(gene_type == "protein_coding")

write.xlsx(AKPvAK, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8A/AKPvAK.xlsx")


AKPSvAKP <- lfcShrink(dds_human_transcriptomic, contrast = c("group", "AKPS", "AKP"), type = "ashr") %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_meta %>% dplyr::select(gene_id, gene_name, gene_type), by = c("ensembl_id" = "gene_id")) %>%
  dplyr::filter(gene_type == "protein_coding")

write.xlsx(AKPvAK, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8B/AKPSvAKP.xlsx")

#######################################################
# 8. Immune cell fractions across predicted classes (AK -> AKP -> AKPS)
########################################################
message("8. Preparing CIBERSORT input for predicted (transcriptomic) groups...")

cibersort_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/CIBERSORT"
dir.create(cibersort_dir, showWarnings = FALSE, recursive = TRUE)

AK_predicted_samples   <- pred_results_15$Sample[pred_results_15$Predicted_Class == "AK"]
AKP_predicted_samples  <- pred_results_15$Sample[pred_results_15$Predicted_Class == "AKP"]
AKPS_predicted_samples <- pred_results_15$Sample[pred_results_15$Predicted_Class == "AKPS"]

write_cibersort_input <- function(tpm_mat, samples, out_path) {
  tpm_mat[, colnames(tpm_mat) %in% samples, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("GeneSymbol") %>%
    write.table(out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# NOTE: assumes `coad_tpm_cibersort` (HGNC-symbol log2(TPM+1) matrix

write_cibersort_input(coad_tpm_cibersort, AK_predicted_samples,
                      file.path(cibersort_dir, "COAD_counts_tpm_AK_predicted.txt"))
write_cibersort_input(coad_tpm_cibersort, AKP_predicted_samples,
                      file.path(cibersort_dir, "COAD_counts_tpm_AKP_predicted.txt"))
write_cibersort_input(coad_tpm_cibersort, AKPS_predicted_samples,
                      file.path(cibersort_dir, "COAD_counts_tpm_AKPS_predicted.txt"))

message("Run the three COAD_counts_tpm_*_predicted.txt files through CIBERSORT externally.")
message("Then save results as CIBERSORT_AK_predicted_immunefraction.csv (etc.) in the same folder, and continue below.")

# ---------------------------------------------------------------
# ---------------------------------------------------------------
read_cibersort_raw <- function(file, group_name) {
  read.csv(file) %>%
    dplyr::rename(Sample = Mixture) %>%
    dplyr::select(-P.value, -Correlation, -RMSE, -`Absolute.score..sig.score.`) %>%
    mutate(Group = group_name)
}

cibersort_AK   <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_AK_predicted_immunefraction.csv"), "AK")
cibersort_AKP  <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_AKP_predicted_immunefraction.csv"), "AKP")
cibersort_AKPS <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_AKPS_predicted_immunefraction.csv"), "AKPS")

cibersort_long <- bind_rows(cibersort_AK, cibersort_AKP, cibersort_AKPS) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "cell_type", values_to = "fraction")


cibersort_regroup <- c(
  "B.cells.naive" = "B cells", "B.cells.memory" = "B cells", "Plasma.cells" = "B cells",
  "NK.cells.resting" = "NK cells", "NK.cells.activated" = "NK cells",
  "T.cells.CD4.naive" = "Helper T cells", "T.cells.CD4.memory.resting" = "Helper T cells",
  "T.cells.CD4.memory.activated" = "Helper T cells", "T.cells.follicular.helper" = "Helper T cells",
  "T.cells.CD8" = "CD8 T cells",
  "T.cells.regulatory..Tregs." = "Regulatory T cells",
  "T.cells.gamma.delta" = "Gamma delta T cells",
  "Monocytes" = "Monocytes",
  "Macrophages.M0" = "Macrophages", "Macrophages.M1" = "Macrophages", "Macrophages.M2" = "Macrophages",
  "Dendritic.cells.resting" = "Dendritic cells", "Dendritic.cells.activated" = "Dendritic cells",
  "Mast.cells.resting" = "Mast cells", "Mast.cells.activated" = "Mast cells",
  "Eosinophils" = "Eosinophils", "Neutrophils" = "Neutrophils"
)

cibersort_panel_map <- c(
  "B cells" = "Lymphoid", "NK cells" = "Lymphoid", "Helper T cells" = "Lymphoid",
  "CD8 T cells" = "Lymphoid", "Regulatory T cells" = "Lymphoid", "Gamma delta T cells" = "Lymphoid",
  "Monocytes" = "Myeloid", "Macrophages" = "Myeloid", "Dendritic cells" = "Myeloid",
  "Mast cells" = "Granulocytes", "Eosinophils" = "Granulocytes", "Neutrophils" = "Granulocytes"
)

plot_data_immune <- cibersort_long %>%
  dplyr::filter(cell_type %in% names(cibersort_regroup)) %>%
  mutate(
    cell_type = cibersort_regroup[cell_type],
    panel = cibersort_panel_map[cell_type],
    Group = factor(Group, levels = c("AK", "AKP", "AKPS"))
  ) %>%
  
  group_by(Sample, Group, panel, cell_type) %>%
  summarise(fraction = sum(fraction, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, panel, cell_type) %>%
  summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

# Order cell types by lineage within the single stacked bar, and
# color by lineage block, so the Lymphoid/Myeloid/Granulocyte
# groupings are still visually distinguishable in one combined figure.
cell_type_order <- c(
  "B cells", "NK cells", "Helper T cells", "CD8 T cells", "Regulatory T cells", "Gamma delta T cells",  # Lymphoid
  "Monocytes", "Macrophages", "Dendritic cells",                                                          # Myeloid
  "Mast cells", "Eosinophils", "Neutrophils"                                                               # Granulocytes
)

lineage_colors <- c(
  "B cells" = "#08519c", "NK cells" = "#3182bd", "Helper T cells" = "#6baed6",
  "CD8 T cells" = "#9ecae1", "Regulatory T cells" = "#c6dbef", "Gamma delta T cells" = "#deebf7",
  "Monocytes" = "#e6550d", "Macrophages" = "#fd8d3c", "Dendritic cells" = "#fdbe85",
  "Mast cells" = "#238b45", "Eosinophils" = "#74c476", "Neutrophils" = "#bae4b3"
)

plot_data_immune <- plot_data_immune %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_order))

combined_immune_plot <- ggplot(plot_data_immune, aes(x = Group, y = fraction, fill = cell_type)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = lineage_colors) +
  theme_minimal() +
  labs(x = "Predicted Class", y = "Mean CIBERSORTx Absolute Score", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "immune_fractions_predicted_combined.pdf"),
       plot = combined_immune_plot, width = 8, height = 8)

write.xlsx(plot_data_immune, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8E-G/plot_data_immune.xlsx")

########################################################
# 9. CMS composition across predicted classes (AK -> AKP -> AKPS)
########################################################
message("9. Building CMS composition barplot (predicted classes)...")

CMS_classifications <- read.csv(file.path(genomic_dir, "CMS_classification_results.csv")) %>%
  dplyr::rename(SampleName = X) %>%
  mutate(SampleName = substr(SampleName, 1, 15))

CMS_classifications_predicted <- CMS_classifications %>%
  left_join(pred_results_15[, c("Sample", "Predicted_Class")], by = c("SampleName" = "Sample")) %>%
  dplyr::filter(Predicted_Class %in% c("AK", "AKP", "AKPS")) %>%
  mutate(
    Predicted_Class = factor(Predicted_Class, levels = c("AK", "AKP", "AKPS")),
    RF = replace_na(RF, "Other")
  )

CMS_predicted_grouped <- CMS_classifications_predicted %>%
  group_by(RF, Predicted_Class) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(Predicted_Class) %>%
  mutate(proportion = Freq / sum(Freq)) %>%
  ungroup()

humanStackedCMS_predicted <- ggplot(CMS_predicted_grouped, aes(x = Predicted_Class, y = proportion, fill = RF)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(
    title = "CMS Subtype Composition Across Predicted Classes",
    x = "Predicted Class (Transcriptomic)",
    y = "Fraction of Samples Belonging to CMS Subtype",
    fill = "CMS Subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

ggsave(file.path(output_dir, "Human_barplot_CMS_predicted.pdf"), humanStackedCMS_predicted, width = 14, height = 8)

write.xlsx(CMS_predicted_grouped, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8C/CMS_predicted_grouped.xlsx")
message("Done.")


####### now rerun CMS plot and immune stacked barplot graph based on the new TUMOR classifier ######### 
pred_results_tumor <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/all_pred_results_TumorResults.csv")

pred_results_tumor_15 <- pred_results_tumor %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE)  # guards against duplicate-key row multiplication in any join below

dup_after_trunc <- sum(duplicated(substr(pred_results_tumor$Sample, 1, 15)))
if (dup_after_trunc > 0) {
  warning(sprintf("%d pred_results samples collapsed to a duplicate 15-char ID after truncation; kept first occurrence only.",
                  dup_after_trunc))
}



# 8B. Immune cell fractions across predicted classes (AK -> AKP -> AKPS)
########################################################
message("8. Preparing CIBERSORT input for predicted (transcriptomic) groups...")



cibersort_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/CIBERSORT"
dir.create(cibersort_dir, showWarnings = FALSE, recursive = TRUE)

Polyp_predicted_samples   <- pred_results_tumor_15$Sample[pred_results_tumor_15$Predicted_Class == "Polyps"]
AKPTu_predicted_samples  <- pred_results_tumor_15$Sample[pred_results_tumor_15$Predicted_Class == "AKP_Tu"]
AKPSTu_predicted_samples <- pred_results_tumor_15$Sample[pred_results_tumor_15$Predicted_Class == "AKPS_Tu"]

write_cibersort_input <- function(tpm_mat, samples, out_path) {
  tpm_mat[, colnames(tpm_mat) %in% samples, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("GeneSymbol") %>%
    write.table(out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# NOTE: assumes `coad_tpm_cibersort` (HGNC-symbol log2(TPM+1) matrix -- same as above
write_cibersort_input(coad_tpm_cibersort, Polyp_predicted_samples ,
                      file.path(cibersort_dir, "COAD_counts_tpm_Polyp_predicted.txt"))
write_cibersort_input(coad_tpm_cibersort, AKPTu_predicted_samples,
                      file.path(cibersort_dir, "COAD_counts_tpm_AKPTu_predicted.txt"))
write_cibersort_input(coad_tpm_cibersort, AKPSTu_predicted_samples,
                      file.path(cibersort_dir, "COAD_counts_tpm_AKPSTu_predicted.txt"))

message("Run the three COAD_counts_tpm_*_predicted.txt files through CIBERSORT externally.")
message("Then save results as CIBERSORT_AK_predicted_immunefraction.csv (etc.) in the same folder, and continue below.")

# ---------------------------------------------------------------
# ---------------------------------------------------------------
read_cibersort_raw <- function(file, group_name) {
  read.csv(file) %>%
    dplyr::rename(Sample = Mixture) %>%
    dplyr::select(-P.value, -Correlation, -RMSE, -`Absolute.score..sig.score.`) %>%
    mutate(Group = group_name)
}

cibersort_Polyp   <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_Polyps_predicted_immunefraction.csv"), "Polyp")
cibersort_AKPTu  <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_AKPTu_predicted_immunefraction.csv"), "AKPTu")
cibersort_AKPSTu <- read_cibersort_raw(file.path(cibersort_dir, "CIBERSORT_AKPSTu_predicted_immunefraction.csv"), "AKPSTu")


cibersort_long <- bind_rows(cibersort_Polyp, cibersort_AKPTu, cibersort_AKPSTu ) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "cell_type", values_to = "fraction")


cibersort_regroup <- c(
  "B.cells.naive" = "B cells", "B.cells.memory" = "B cells", "Plasma.cells" = "B cells",
  "NK.cells.resting" = "NK cells", "NK.cells.activated" = "NK cells",
  "T.cells.CD4.naive" = "Helper T cells", "T.cells.CD4.memory.resting" = "Helper T cells",
  "T.cells.CD4.memory.activated" = "Helper T cells", "T.cells.follicular.helper" = "Helper T cells",
  "T.cells.CD8" = "CD8 T cells",
  "T.cells.regulatory..Tregs." = "Regulatory T cells",
  "T.cells.gamma.delta" = "Gamma delta T cells",
  "Monocytes" = "Monocytes",
  "Macrophages.M0" = "Macrophages", "Macrophages.M1" = "Macrophages", "Macrophages.M2" = "Macrophages",
  "Dendritic.cells.resting" = "Dendritic cells", "Dendritic.cells.activated" = "Dendritic cells",
  "Mast.cells.resting" = "Mast cells", "Mast.cells.activated" = "Mast cells",
  "Eosinophils" = "Eosinophils", "Neutrophils" = "Neutrophils"
)

cibersort_panel_map <- c(
  "B cells" = "Lymphoid", "NK cells" = "Lymphoid", "Helper T cells" = "Lymphoid",
  "CD8 T cells" = "Lymphoid", "Regulatory T cells" = "Lymphoid", "Gamma delta T cells" = "Lymphoid",
  "Monocytes" = "Myeloid", "Macrophages" = "Myeloid", "Dendritic cells" = "Myeloid",
  "Mast cells" = "Granulocytes", "Eosinophils" = "Granulocytes", "Neutrophils" = "Granulocytes"
)

plot_data_immune <- cibersort_long %>%
  dplyr::filter(cell_type %in% names(cibersort_regroup)) %>%
  mutate(
    cell_type = cibersort_regroup[cell_type],
    panel = cibersort_panel_map[cell_type],
    Group = factor(Group, levels = c("Polyp", "AKPTu", "AKPSTu"))
  ) %>%
  
  group_by(Sample, Group, panel, cell_type) %>%
  summarise(fraction = sum(fraction, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, panel, cell_type) %>%
  summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

# Order cell types by lineage within the single stacked bar, and
# color by lineage block, so the Lymphoid/Myeloid/Granulocyte
# groupings are still visually distinguishable in one combined figure.
cell_type_order <- c(
  "B cells", "NK cells", "Helper T cells", "CD8 T cells", "Regulatory T cells", "Gamma delta T cells",  # Lymphoid
  "Monocytes", "Macrophages", "Dendritic cells",                                                          # Myeloid
  "Mast cells", "Eosinophils", "Neutrophils"                                                               # Granulocytes
)

lineage_colors <- c(
  "B cells" = "#08519c", "NK cells" = "#3182bd", "Helper T cells" = "#6baed6",
  "CD8 T cells" = "#9ecae1", "Regulatory T cells" = "#c6dbef", "Gamma delta T cells" = "#deebf7",
  "Monocytes" = "#e6550d", "Macrophages" = "#fd8d3c", "Dendritic cells" = "#fdbe85",
  "Mast cells" = "#238b45", "Eosinophils" = "#74c476", "Neutrophils" = "#bae4b3"
)

plot_data_immune <- plot_data_immune %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_order))

combined_immune_plot <- ggplot(plot_data_immune, aes(x = Group, y = fraction, fill = cell_type)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = lineage_colors) +
  theme_minimal() +
  labs(x = "Predicted Class", y = "Mean CIBERSORTx Absolute Score", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "immune_fractions_predictedTu_combined.pdf"),
       plot = combined_immune_plot, width = 8, height = 8)

write.xlsx(plot_data_immune, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8B/plot_data_immune_tu.xlsx")

########################################################
# 9B. CMS composition across predicted classes (Polyps -> AKP_Tu -> AKPS_Tu)
########################################################
message("9. Building CMS composition barplot (predicted classes)...")

CMS_classifications <- read.csv(file.path(genomic_dir, "CMS_classification_results.csv")) %>%
  dplyr::rename(SampleName = X) %>%
  mutate(SampleName = substr(SampleName, 1, 15))

CMS_classifications_predicted <- CMS_classifications %>%
  left_join(pred_results_tumor_15[, c("Sample", "Predicted_Class")], by = c("SampleName" = "Sample")) %>%
  dplyr::filter(Predicted_Class %in% c("Polyps", "AKP_Tu", "AKPS_Tu")) %>%
  mutate(
    Predicted_Class = factor(Predicted_Class, levels = c("Polyps", "AKP_Tu", "AKPS_Tu")),
    RF = replace_na(RF, "Other")
  )

CMS_predicted_grouped <- CMS_classifications_predicted %>%
  group_by(RF, Predicted_Class) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(Predicted_Class) %>%
  mutate(proportion = Freq / sum(Freq)) %>%
  ungroup()

humanStackedCMS_predicted <- ggplot(CMS_predicted_grouped, aes(x = Predicted_Class, y = proportion, fill = RF)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(
    title = "CMS Subtype Composition Across Predicted Classes",
    x = "Predicted Class (Transcriptomic)",
    y = "Fraction of Samples Belonging to CMS Subtype",
    fill = "CMS Subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

ggsave(file.path(output_dir, "Human_barplot_CMS_predicted_Tu.pdf"), humanStackedCMS_predicted, width = 14, height = 8)

write.xlsx(CMS_predicted_grouped, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8D/CMS_predicted_grouped_tu.xlsx")
message("Done.")


####################### Sup Fig 9 #################### 
output_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures"

# ---------------------------------------------------------------
--------------------
pred_results_immune <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/pred_results_IMMUNE_NAIVE_ORGANOIDS.csv") %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  dplyr::select(Sample, Immune_Class = Predicted_Class)

pred_results_tumor <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/pred_results_IMMUNE_EXPOSED_TUMORS.csv") %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  dplyr::select(Sample, Tumor_Class = Predicted_Class)

class_compare <- inner_join(pred_results_immune, pred_results_tumor, by = "Sample")

n_immune_only <- sum(!pred_results_immune$Sample %in% pred_results_tumor$Sample)
n_tumor_only  <- sum(!pred_results_tumor$Sample %in% pred_results_immune$Sample)
if (n_immune_only > 0 || n_tumor_only > 0) {
  message(sprintf("%d samples only in immune-naive predictions, %d samples only in tumor predictions -- excluded from the confusion matrix (only samples present in BOTH are compared).",
                  n_immune_only, n_tumor_only))
}

# ---------------------------------------------------------------
# 2. Build the confusion matrix.
# ---------------------------------------------------------------
order_with_none_last <- function(labels) {
  labels <- sort(unique(labels))
  if ("None" %in% labels) labels <- c(setdiff(labels, "None"), "None")
  labels
}

immune_labels <- order_with_none_last(class_compare$Immune_Class)
tumor_labels  <- order_with_none_last(class_compare$Tumor_Class)

table_compare <- table(
  factor(class_compare$Immune_Class, levels = immune_labels),
  factor(class_compare$Tumor_Class, levels = tumor_labels)
)

# ---------------------------------------------------------------
# 3. Plot
# ---------------------------------------------------------------
png(file.path(output_dir, "ConfusionMatrix_ImmuneNaive_vs_Tumor.png"), width = 6, height = 6, units = "in", res = 300)
ph <- pheatmap::pheatmap(
  table_compare, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
  color = colorRampPalette(c("white", "steelblue"))(50), fontsize = 15, fontsize_number = 13,
  angle_col = 270, number_format = "%.0f", cellwidth = 50, cellheight = 50, silent = TRUE
)
grid.newpage()
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.draw(ph$gtable)
grid.text("Tumor Classifier Prediction", x = 0.5, y = unit(0.02, "npc"), gp = gpar(fontsize = 12))
grid.text("Immune-Naive Organoid Classifier Prediction", x = unit(0.02, "npc"), y = 0.5, rot = 90, gp = gpar(fontsize = 12))
dev.off()

message(sprintf("Confusion matrix built from %d samples with predictions from both classifiers.", nrow(class_compare)))
print(table_compare)

write.xlsx(table_compare, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S7/tableCompare.xlsx")


###################### now add the unclassified CMS analysis to respond to reviewer 3 ############## 
output_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output"
genomic_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data"

raw_counts_unique <- assay(coad_data)[, sample_info_unique$full_id]
colnames(raw_counts_unique) <- sample_info_unique$base_id

CMS_classifications <- read.csv(file.path(genomic_dir, "CMS_classification_results.csv")) %>%
  dplyr::rename(SampleName = X) %>%
  mutate(SampleName = substr(SampleName, 1, 15))

#### make CMS metadata now 

cms_sample_metadata <- data.frame(sample_id = colnames(raw_counts_unique), stringsAsFactors = FALSE) %>%
  inner_join(CMS_classifications, by = c("sample_id" = "SampleName")) %>%  # inner_join: only samples the classifier actually scored
  mutate(CMS_status = if_else(is.na(RF), "Unclassified", "Classified"))

n_never_scored <- sum(!colnames(raw_counts_unique) %in% CMS_classifications$SampleName)
message(sprintf("%d raw-count samples were never scored by the CMS classifier and are excluded from this comparison entirely (not counted as Unclassified).",
                n_never_scored))
message(sprintf("Of %d scored samples: %d Classified, %d Unclassified.",
                nrow(cms_sample_metadata),
                sum(cms_sample_metadata$CMS_status == "Classified"),
                sum(cms_sample_metadata$CMS_status == "Unclassified")))


message("Building DGE (Unclassified vs. Classified)...")

dds_cms_unclassified <- DESeqDataSetFromMatrix(
  countData = raw_counts_unique[, cms_sample_metadata$sample_id],
  colData   = cms_sample_metadata,
  design    = ~ CMS_status
)
stopifnot(identical(colnames(dds_cms_unclassified), cms_sample_metadata$sample_id))
dds_cms_unclassified <- DESeq(dds_cms_unclassified)

res_cms <- results(dds_cms_unclassified, contrast = c("CMS_status", "Unclassified", "Classified")) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_meta %>% dplyr::select(gene_id, gene_name), by = c("ensembl_id" = "gene_id")) %>%
  dplyr::filter(!is.na(stat), !is.na(gene_name))

write.csv(res_cms, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S8/res_cms.csv", row.names = TRUE)

volcano_unclassified_vs_classified <- make_volcano_human(
  dds_cms_unclassified,
  c("CMS_status", "Unclassified", "Classified"),
  "CMS-Unclassified vs. CMS-Classified",
  show_legend = TRUE
)

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures/unclassifiedVClassified.pdf",
       volcano_unclassified_vs_classified, width = 10, height = 8)

############### do GSEA now ############
res_ranked <- res_cms %>%
  arrange(desc(stat)) %>%
  dplyr::select(gene_name, stat)

rnk_path <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S8/unClassvClass.rnk"
write.table(res_ranked, file = rnk_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

############### make heatmap ########## 
run_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/UnClassvClass.GseaPreranked.1789325916822"

neg_file <- list.files(run_dir, pattern = "^gsea_report_for_na_neg_.*\\.tsv$", full.names = TRUE)
pos_file <- list.files(run_dir, pattern = "^gsea_report_for_na_pos_.*\\.tsv$", full.names = TRUE)
stopifnot(length(neg_file) == 1, length(pos_file) == 1)

gsea_unclass_vs_class <- bind_rows(
  read_and_clean(neg_file),
  read_and_clean(pos_file)
) %>%
  dplyr::rename(pathway = NAME, padj = `FDR q-val`) %>%
  filter(!is.na(NES))

gsea_results_list <- list("UnClassvClass" = gsea_unclass_vs_class)
comparison_labels <- c("UnClassvClass" = "Unclassified vs. Classified")

heatmap_immune <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list, comparison_labels, immune_pathways,
  xlsx_path = file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S8", "UnClassvClass_immuneNES.xlsx"),
  save_path = file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures", "hallMarkHeatMap_UnClassvClass_immune.pdf"),
  save_width = 14, save_height = 10, pval_cutoff = 0.05
)

heatmap_oncogenic <- plot_full_hallmark_heatmap_with_priority(
  gsea_results_list, comparison_labels, oncogenic_pathways,
  xlsx_path = file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S8", "UnClassvClass_oncogenicNES.xlsx"),
  save_path = file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures", "hallMarkHeatMap_UnClassvClass_oncogenic.pdf"),
  save_width = 14, save_height = 10, pval_cutoff = 0.05,
  strip_text_size = 9
)

#######now look at A-AK-AKP-AKPS with respect to this unclassified group ########

# ---------------------------------------------------------------
gene_cols <- c("APC", "KRAS", "TP53", "SMAD4")

# ---------------------------------------------------------------
# 1. Load mutation calls
# ---------------------------------------------------------------
tcga_mutation <- read_delim(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 1/Figure A/upsetTCGA.txt",
  delim = "\t"
) %>%
  mutate(
    SampleName = substr(`Sample ID`, 1, 15),
    across(all_of(gene_cols), ~ as.integer(. == "Yes"))
  )

# ---------------------------------------------------------------
# 2. Mutation-combination group (A/AK/AKP/AKPS) 
# ---------------------------------------------------------------
tcga_mutation <- tcga_mutation %>%
  rowwise() %>%
  mutate(mutation_group_raw = paste0(
    if_else(APC == 1, "A", ""),
    if_else(KRAS == 1, "K", ""),
    if_else(TP53 == 1, "P", ""),
    if_else(SMAD4 == 1, "S", "")
  )) %>%
  ungroup() %>%
  mutate(mutation_group = if_else(mutation_group_raw %in% c("A", "AK", "AKP", "AKPS"),
                                  mutation_group_raw, "Other"))

mutation_cms <- cms_sample_metadata %>%
  dplyr::select(sample_id, CMS_status) %>%
  inner_join(
    tcga_mutation %>% dplyr::select(SampleName, mutation_group),
    by = c("sample_id" = "SampleName")
  )

n_no_mutation_data <- sum(!cms_sample_metadata$sample_id %in% tcga_mutation$SampleName)
message(sprintf("%d CMS-scored samples had no matching row in the mutation file and are excluded from this comparison.",
                n_no_mutation_data))
message(sprintf("Comparing mutation group composition across %d samples: %d Classified, %d Unclassified.",
                nrow(mutation_cms),
                sum(mutation_cms$CMS_status == "Classified"),
                sum(mutation_cms$CMS_status == "Unclassified")))

# ---------------------------------------------------------------
# 3. Stacked bar plot: mutation group composition by CMS status
# ---------------------------------------------------------------
mutation_group_summary <- mutation_cms %>%
  mutate(mutation_group = factor(mutation_group, levels = c("A", "AK", "AKP", "AKPS", "Other"))) %>%
  dplyr::count(CMS_status, mutation_group) %>%
  group_by(CMS_status) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

mutation_group_barplot <- ggplot(mutation_group_summary, aes(x = CMS_status, y = proportion, fill = mutation_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(
    title = "Mutation Group Composition: CMS-Classified vs. Unclassified",
    x = "CMS Status", y = "Proportion of Samples", fill = "Mutation Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

ggsave(file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures/", "CMS_unclassified_mutation_group_barplot.pdf"),
       mutation_group_barplot, width = 8, height = 8)
write.csv(mutation_group_summary, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S8", "CMS_unclassified_mutation_group_summary.csv"), row.names = FALSE)

########## now examine CISR groups per ?reviwer 2 suggestion #########
cris_input <- as.matrix(coad_tpm_cibersort)
# 1. Run CRIS classification via CMScaller's nearest-template-
#    prediction (ntp()) using templates.CRIS.

message("Running CRIS classification (nearest template prediction)...")
templates_cris_symbol <- CMScaller::templates.CRIS
templates_cris_symbol$probe <- templates_cris_symbol$symbol

cris_result <- ntp(emat = cris_input, templates = templates_cris_symbol, doPlot = FALSE)

cris_classifications <- cris_result %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleName") %>%
  mutate(
    SampleName = substr(SampleName, 1, 15),
    CRIS_subtype = if_else(is.na(prediction), "Unclassified", as.character(prediction))
  )

write.csv(cris_classifications, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S9", "CRIS_classification_results.csv"), row.names = FALSE)

pred_immune_naive <- read.csv(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 7/pred_results_IMMUNE_NAIVE_ORGANOIDS.csv"
) %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE)

pred_immune_exposed <- read.csv(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/pred_results_IMMUNE_EXPOSED_TUMORS.csv"
) %>%
  mutate(Sample = substr(Sample, 1, 15)) %>%
  distinct(Sample, .keep_all = TRUE)


#### resuable function 
build_cris_composition_plot <- function(pred_df, class_levels, title, save_name) {
  merged <- pred_df %>%
    dplyr::select(Sample, Predicted_Class) %>%
    inner_join(cris_classifications %>% dplyr::select(SampleName, CRIS_subtype), by = c("Sample" = "SampleName")) %>%
    filter(Predicted_Class %in% class_levels) %>%
    mutate(Predicted_Class = factor(Predicted_Class, levels = class_levels))
  
  n_unmatched <- sum(!pred_df$Sample %in% cris_classifications$SampleName)
  message(sprintf("  %s: %d/%d samples matched to a CRIS classification (%d unmatched).",
                  title, nrow(merged), nrow(pred_df), n_unmatched))
  
  summary_df <- merged %>%
    dplyr::count(Predicted_Class, CRIS_subtype) %>%
    group_by(Predicted_Class) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  p <- ggplot(summary_df, aes(x = Predicted_Class, y = proportion, fill = CRIS_subtype)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    labs(title = title, x = "Predicted Class", y = "Proportion of Samples", fill = "CRIS Subtype") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())
  
  ggsave(file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures", save_name), p, width = 8, height = 8)
  list(plot = p, data = summary_df)
}

message("Building CRIS composition barplots...")

result_immune_naive <- build_cris_composition_plot(
  pred_immune_naive, c("AK", "AKP", "AKPS"),
  "CRIS Subtype Composition: AK -> AKP -> AKPS",
  "CRIS_composition_AK_AKP_AKPS.pdf"
)

result_immune_exposed <- build_cris_composition_plot(
  pred_immune_exposed, c("Polyps", "AKP_Tu", "AKPS_Tu"),
  "CRIS Subtype Composition: Polyps -> AKP_Tu -> AKPS_Tu",
  "CRIS_composition_Polyps_AKPTu_AKPSTu.pdf"
)

message("Done.")
