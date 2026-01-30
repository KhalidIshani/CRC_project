suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggalluvial)
  library(tidyr)
  library(purrr)
})

# ================================================================
# --- PARAMETERS (EDIT ONLY THIS BLOCK) ---
# ================================================================
DATA_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/genomic"
OUTPUT_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/genomic"

# ✅ You can skip any of these by setting to NULL or ""
A_file    <- NULL
AK_file   <- file.path(DATA_DIR, "AK_genomic.csv")
AKP_file  <- file.path(DATA_DIR, "AKP_genomic.csv")
AKPS_file <- file.path(DATA_DIR, "AKPS_genomic.csv")

# ================================================================
# --- 1. Download and Prepare TCGA Data ---
# ================================================================
download_tcga_data <- function() {
  message("Downloading TCGA-COAD data...")
  query <- GDCquery(project = 'TCGA-COAD',
                    data.category = 'Transcriptome Profiling',
                    data.type = 'Gene Expression Quantification',
                    workflow.type = 'STAR - Counts')
  GDCdownload(query)
  GDCprepare(query, summarizedExperiment = TRUE)
}

# ================================================================
# --- 2. TPM Conversion ---
# ================================================================
convert_to_tpm <- function(se_obj) {
  message("Converting counts to log2(TPM + 1)...")
  
  meta <- as.data.frame(se_obj@rowRanges) %>%
    dplyr::select(gene_id, width) %>%
    mutate(width_kb = width / 1000)
  
  counts <- as.data.frame(assay(se_obj)) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(meta, by = c("ensembl_id" = "gene_id")) %>%
    filter(!is.na(width_kb) & width_kb > 0)
  
  counts_mat <- as.matrix(counts[, !(names(counts) %in% c("ensembl_id", "width", "width_kb"))])
  rownames(counts_mat) <- counts$ensembl_id
  
  counts_per_kb <- sweep(counts_mat, 1, counts$width_kb, "/")
  tpm <- sweep(counts_per_kb, 2, colSums(counts_per_kb), "/") * 1e6
  log2_tpm <- log2(tpm + 1)
  
  list(tpm = tpm, log2_tpm = log2_tpm)
}

# ================================================================
# --- 3. Filter and Rename Samples ---
# ================================================================
filter_samples <- function(tpm_matrix, sample_lists) {
  message("Filtering and renaming samples...")
  
  subsets <- list()
  for (grp in names(sample_lists)) {
    samples <- sample_lists[[grp]]
    
    # Skip if samples vector is empty or NULL
    if (is.null(samples) || length(samples) == 0) {
      warning("⚠️  Skipping ", grp, " — no samples provided.")
      next
    }
    
    # Match sample IDs (trimmed to first 15 characters for TCGA consistency)
    tpm_cols <- substr(colnames(tpm_matrix), 1, 15)
    valid_samples <- colnames(tpm_matrix)[tpm_cols %in% substr(samples, 1, 15)]
    
    if (length(valid_samples) == 0) {
      warning("⚠️  Skipping ", grp, " — no matching samples found in TPM matrix.")
      next
    }
    
    subsets[[grp]] <- tpm_matrix[, valid_samples, drop = FALSE]
    message("✅ ", grp, ": ", length(valid_samples), " samples matched.")
  }
  
  return(subsets)
}

# ================================================================
# --- 4. Prepare CIBERSORT Input ---
# ================================================================
# FAST + SAFE version of prepare_cibersort_input()
prepare_cibersort_input <- function(tpm_matrix, sample_lists, outdir) {
  message("\n🚀 Preparing CIBERSORT input files...")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(data.table)
    library(tibble)
  })
  
  ## --- Step 1: Clean Ensembl IDs ---
  message("Step 1: Cleaning Ensembl IDs...")
  ensembl_clean <- sub("\\..*$", "", rownames(tpm_matrix))
  
  ## --- Step 2: Fast Ensembl→HGNC mapping ---
  message("Step 2: Mapping Ensembl → HGNC symbols (vectorized)...")
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ensembl_clean),
    columns = c("ENSEMBL", "SYMBOL"),
    keytype = "ENSEMBL"
  )
  hgnc <- mapping$SYMBOL[match(ensembl_clean, mapping$ENSEMBL)]
  valid <- !is.na(hgnc)
  tpm_hgnc <- tpm_matrix[valid, , drop = FALSE]
  hgnc <- hgnc[valid]
  
  ## --- Step 3: Average TPM by gene symbol ---
  message("Step 3: Averaging expression by gene symbol (fast data.table)...")
  dt <- as.data.table(tpm_hgnc)
  dt[, HGNC_Symbol := hgnc]
  tpm_avg <- dt[, lapply(.SD, mean, na.rm = TRUE), by = HGNC_Symbol]
  
  # Reformat for downstream use
  tpm_avg_mat <- as.data.frame(tpm_avg)
  rownames(tpm_avg_mat) <- tpm_avg_mat$HGNC_Symbol
  tpm_avg_mat$HGNC_Symbol <- NULL
  
  ## --- Step 4: Write per-group CIBERSORT files ---
  message("Step 4: Writing group files for CIBERSORT...")
  for (grp in names(sample_lists)) {
    samples <- colnames(sample_lists[[grp]])
    valid_samples <- intersect(samples, colnames(tpm_avg_mat))
    
    if (length(valid_samples) == 0) {
      message("  ⚠️  Skipping ", grp, " — no matching samples found.")
      next
    }
    
    # Optionally skip single-sample groups (CIBERSORTx may fail on n=1)
    if (length(valid_samples) < 2) {
      message("  ⚠️  Skipping ", grp, " — only one valid sample.")
      next
    }
    
    sub <- tpm_avg_mat[, valid_samples, drop = FALSE]
    sub_df <- tibble::rownames_to_column(as.data.frame(sub), "GeneSymbol")
    
    out_file <- file.path(outdir, paste0("COAD_counts_tpm_", grp, ".txt"))
    write.table(sub_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("  ✅ Wrote ", grp, " (", ncol(sub_df) - 1, " samples) → ", out_file)
  }
  
  message("\n🎉 Done! CIBERSORT input files are in: ", normalizePath(outdir))
}

# ================================================================
# --- 5. Wilcoxon Tests for CIBERSORT Results ---
# ================================================================
compare_cibersort_groups <- function(file1, file2, label1, label2, out_file) {
  if (!file.exists(file1) || !file.exists(file2)) {
    message("Skipping comparison ", label1, " vs ", label2, " (missing file)")
    return(invisible(NULL))
  }
  
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  cols <- setdiff(names(df1), c("Mixture", "P.value", "Correlation", "RMSE", "Absolute.score..sig.score."))
  
  res <- lapply(cols, function(cell_type) {
    g1 <- df1[[cell_type]]
    g2 <- df2[[cell_type]]
    test <- tryCatch(wilcox.test(g1, g2, paired = FALSE, alternative = "less"),
                     error = function(e) list(p.value = NA))
    c(mean1 = mean(g1, na.rm = TRUE),
      mean2 = mean(g2, na.rm = TRUE),
      p = test$p.value)
  })
  
  df_res <- as.data.frame(do.call(rbind, res))
  df_res$cell_type <- cols
  names(df_res)[1:2] <- c(paste0("mean_", label1), paste0("mean_", label2))
  df_res$FDR <- p.adjust(df_res$p, "BH")
  write.csv(df_res, out_file, row.names = FALSE)
  message("Written: ", out_file)
}

# ================================================================
# --- 6. Alluvial Plot ---
# ================================================================
plot_alluvial <- function(cibersort_files, groups, outfile) {
  existing <- cibersort_files[file.exists(cibersort_files)]
  groups <- groups[file.exists(cibersort_files)]
  if (length(existing) < 2) {
    message("Not enough groups for alluvial plot.")
    return(invisible(NULL))
  }
  
  read_cibersort <- function(file, group) {
    df <- read.csv(file, check.names = FALSE)
    
    # Columns to remove (if they exist)
    cols_to_remove <- c(
      "Mixture",
      "P.value",
      "P-value",
      "Correlation",
      "RMSE",
      "Absolute.score..sig.score.",
      "Absolute.score",
      "Absolute score (sig.score)"
    )
    
    df <- df %>%
      dplyr::select(!any_of(cols_to_remove))
    
    df$Group <- group
    return(df)
  }
  
  df <- purrr::map2_dfr(existing, groups, read_cibersort)
  
  df <- df %>%
    transmute(
      B_cells       = `B cells naive` + `B cells memory`,
      Plasma_cells  = `Plasma cells`,
      T_cells       = `T cells CD8` + `T cells CD4 naive` +
        `T cells CD4 memory resting` + `T cells CD4 memory activated` +
        `T cells follicular helper` + `T cells regulatory (Tregs)` +
        `T cells gamma delta`,
      NK_cells      = `NK cells resting` + `NK cells activated`,
      Monocytes     = `Monocytes`,
      Macrophages   = `Macrophages M0` + `Macrophages M1` + `Macrophages M2`,
      Dendritic_cells = `Dendritic cells resting` + `Dendritic cells activated`,
      Mast_cells      = `Mast cells resting` + `Mast cells activated`,
      Eosinophils     = `Eosinophils`,
      Neutrophils     = `Neutrophils`,
      Group = `Group`
    )
  
  df_long <- df %>%
    pivot_longer(-Group, names_to = "CellType", values_to = "Fraction") %>%
    group_by(Group, CellType) %>%
    summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")
  
  df_long$Group <- factor(df_long$Group, levels = groups)
  
  p <- ggplot(df_long,
              aes(x = Group, stratum = CellType, alluvium = CellType,
                  y = Fraction, fill = CellType)) +
    geom_flow(stat = "alluvium", alpha = 0.7) +
    geom_stratum(width = 0.25) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Immune Cell Counts Across Tumor Progression",
      y = "Immune Cell Count",       # changed here
      x = "Group"
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 18),  # increase legend font
      axis.text.x = element_text(size = 14),  # increase x-axis tick font
      axis.title.x = element_text(size = 18), # x-axis title font
      axis.title.y = element_text(size = 16)  # y-axis title font
    )
  
  ggsave(outfile, p, width = 14, height = 8)
  message("Alluvial plot saved to: ", outfile)
}

# ================================================================
# --- 7. Run the Pipeline ---
# ================================================================
COAD_data_obj <- download_tcga_data()
tpm_results <- convert_to_tpm(COAD_data_obj)

sample_files <- list(
  AK = AK_file,
  AKP = AKP_file,
  AKPS = AKPS_file
)

sample_lists <- lapply(sample_files, function(path) {
  readr::read_csv(path, show_col_types = FALSE)[[3]]  # assumes IDs in first column
})
names(sample_lists) <- names(sample_files)

colnames(tpm_results$tpm) <- substr(colnames(tpm_results$tpm), 1, 15) 
colnames(tpm_results$log2_tpm) <- substr(colnames(tpm_results$log2_tpm), 1, 15) 


sample_subsets <- filter_samples(tpm_results$log2_tpm, sample_lists)


prepare_cibersort_input(tpm_results$tpm, sample_subsets, file.path(DATA_DIR, "CIBERSORT"))


DATA_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/genomic"

# After external CIBERSORT run
groups_present <- names(sample_subsets)
if (length(groups_present) > 1) {
  for (i in seq_len(length(groups_present) - 1)) {
    g1 <- groups_present[i]
    g2 <- groups_present[i + 1]
    compare_cibersort_groups(
      file.path(DATA_DIR, paste0("CIBERSORT_", g1, "_immunefraction.csv")),
      file.path(DATA_DIR, paste0("CIBERSORT_", g2, "_immunefraction.csv")),
      g1, g2,
      file.path(DATA_DIR, paste0("CIBERSORT_results_", g1, "_vs_", g2, ".csv"))
    )
  }
}

plot_alluvial(
  cibersort_files = file.path(DATA_DIR,
                              paste0("CIBERSORT_", groups_present, "_immunefraction.csv")),
  groups = groups_present,
  outfile = file.path(OUTPUT_DIR, "Human_alluvialplot_transcriptomic_substypes.pdf")
)


##### remake genomic alluvial plot but only with AK, AKP, AKPS

DATA_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Genomic"
OUTPUT_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Genomic"
groups_present <- c("AK","AKP","AKPS")

plot_alluvial(
  cibersort_files = file.path(DATA_DIR,
                              paste0("CIBERSORT_", groups_present, "_immunefraction.csv")),
  groups = groups_present,
  outfile = file.path(OUTPUT_DIR, "Human_alluvialplot_genomic_Dec13.pdf")
)



##### remake transcriptomic alluvial plot but only with AK, AKP, AKPS
DATA_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Transcriptomic"
OUTPUT_DIR <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Transcriptomic"
groups_present <- c("AK","AKP","AKPS")


plot_alluvial(
  cibersort_files = file.path(DATA_DIR,
                              paste0("CIBERSORT_", groups_present, "_immunefraction.csv")),
  groups = groups_present,
  outfile = file.path(OUTPUT_DIR, "Human_alluvialplot_transcriptomic_Dec12.pdf")
)

