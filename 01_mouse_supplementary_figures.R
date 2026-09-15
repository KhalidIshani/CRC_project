#===============================================================
# MOUSE COLON CANCER RNA-SEQ -- SUPPLEMENTARY FIGURES
# Includes AK_Tu (different sequencing batch) and the n=1
# metastatic-site groups (AKPS_Met, AKPS_TdLN, AKPS_LungMet).
# Run this script BEFORE 02_mouse_main_figures.R -- the main
# figures script reuses some objects 
#===============================================================

library(tidyverse)
library(DESeq2)
library(umap)
library(UpSetR)
library(RColorBrewer)
library(corrplot)
library(pheatmap)
library(reshape2)
library(DGEobj.utils)
library(janitor)
library(biomaRt)
library(ggalluvial)
library(EnhancedVolcano)
library(gridExtra)
library(cowplot)
library(immunedeconv)
library(data.table)
library(ComplexUpset)
library(dynamicTreeCut)
library(WGCNA)
library(nnet)
library(ggpubr)
library(AnnotationHub)
library(ensembldb)
library(CMScaller)
library(tibble)
library(ggplot2)
if (!requireNamespace("AnnotationHub", quietly = TRUE) || !requireNamespace("ensembldb", quietly = TRUE)) {
  BiocManager::install(c("AnnotationHub", "ensembldb"), update = FALSE, ask = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("icbi-lab/immunedeconv")
}
if (!requireNamespace("MmCMS", quietly = TRUE)) {
  devtools::install_github("MolecularPathologyLab/MmCMS")
}
library(MmCMS)

# ---------------------------------------------------------------
# Config
# ---------------------------------------------------------------
project_dir      <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project"
data_output_dir  <- file.path(project_dir, "Output")
supp_fig_dir     <- file.path(project_dir, "Genome Medicine Revision")
rnk_dir          <- file.path(project_dir, "GSEA/Preranked Genesets")
gsea_out_dir     <- file.path(project_dir, "GSEA/Output")
gsea_res_dir     <- file.path(project_dir, "GSEA/Results")

for (d in c(data_output_dir, supp_fig_dir, rnk_dir, gsea_res_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

counts_file_1 <- "/Users/khalidishani/Downloads/subread_counts.txt"
counts_file_2 <- "/Users/khalidishani/Downloads/subread_counts-2.txt"
metadata_file <- file.path(project_dir, "all_samples_meta.csv")

wanted_samples <- c(
  "AK_Tu_1", "AK_Tu_2", "AK_Tu_3", "NormalColon_002", "NormalColon_003", "NormalColon_004",
  "Polyps_2", "Polyps_3", "Polyp_1",
  "AKPS_LungMet-1_Rep_1", "AKPS_LungMet-1_Rep_2", "AKPS_LungMet-1_Rep_3",
  "AKPS_Met-1_Rep_1", "AKPS_Met-1_Rep_2", "AKPS_Met-1_Rep_3",
  "AKPS_TdLN-1_Rep_1", "AKPS_TdLN-1_Rep_2", "AKPS_TdLN-1_Rep_3",
  "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3",
  "AKP_Tu-1", "AKP_Tu-2_Rep_1", "AKP_Tu-2_Rep_2",
  "AK_2", "AKPS_1", "AK_3", "AKPS_005", "AKP_005", "AKP_001", "AK_1", "AKPS_006", "AKP_004",
  "AKPS_TuOrganoid-2_Rep_2", "AKPS_TuOrganoid-2_Rep_3", "AKP_TuOrganoid-1",
  "AKP_TuOrganoid-2_Rep_1", "AKPS_TuOrganoid-2_Rep_1", "AKP_TuOrganoid-2_Rep_2"
)



# ---------------------------------------------------------------
# 1. Load and merge the two subread count files
#    (join on Geneid, not bind_cols)
# ---------------------------------------------------------------
message("1. Loading and merging raw count files...")

exp1 <- as.data.frame(fread(counts_file_1, sep = "\t"))
exp2 <- as.data.frame(fread(counts_file_2, sep = "\t"))
stopifnot(all(c("Geneid", "Length") %in% names(exp1)))
stopifnot(all(c("Geneid", "Length") %in% names(exp2)))

total_exp <- full_join(exp1, exp2, by = "Geneid", suffix = c(".f1", ".f2"))

length_mismatch <- with(total_exp, !is.na(Length.f1) & !is.na(Length.f2) & Length.f1 != Length.f2)
if (any(length_mismatch, na.rm = TRUE)) {
  warning(sum(length_mismatch), " genes have mismatched Length between the two count files.")
}
total_exp <- total_exp %>%
  mutate(Length = coalesce(Length.f1, Length.f2)) %>%
  dplyr::select(-Length.f1, -Length.f2)

missing_samples <- setdiff(wanted_samples, names(total_exp))
if (length(missing_samples) > 0) {
  stop("These requested samples are not present in either count file: ",
       paste(missing_samples, collapse = ", "))
}

# ---------------------------------------------------------------
# 2. Subset to samples of interest, drop all-zero genes
# ---------------------------------------------------------------
message("2. Subsetting to samples of interest...")

mouse_counts_ensembl <- total_exp[, c(wanted_samples, "Geneid", "Length")]
all_zero <- rowSums(mouse_counts_ensembl[, wanted_samples] == 0) == length(wanted_samples)
mouse_counts_ensembl <- mouse_counts_ensembl[!all_zero, ]
message(sprintf("   Dropped %d genes with zero counts across all samples (%d remain).",
                sum(all_zero), nrow(mouse_counts_ensembl)))


# ---------------------------------------------------------------
# 3. Map Ensembl gene IDs to MGI gene symbols

# ---------------------------------------------------------------
message("3. Mapping Ensembl IDs to MGI symbols...")

# AnnotationHub/ensembldb
if (!requireNamespace("AnnotationHub", quietly = TRUE) || !requireNamespace("ensembldb", quietly = TRUE)) {
  BiocManager::install(c("AnnotationHub", "ensembldb"), update = FALSE, ask = FALSE)
}
library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()
mouse_ensdb_records <- query(ah, pattern = c("Mus musculus", "EnsDb", "110"))


mouse_ensdb_meta <- as.data.frame(mcols(mouse_ensdb_records))
print(mouse_ensdb_meta[, c("title", "species", "genome")])

# Standard reference assembly only (excludes MGP strain-specific
# records, which show a strain name in genome/title instead of
# GRCm38/GRCm39)

standard_idx <- which(mouse_ensdb_meta$species == "Mus musculus" &
                        grepl("^GRCm3[89]$", mouse_ensdb_meta$genome))
stopifnot(length(standard_idx) >= 1)
if (length(standard_idx) > 1) {
  message("Multiple standard-reference records matched; using the first. Check mouse_ensdb_meta[standard_idx, ] if unsure.")
}
mouse_ensdb <- mouse_ensdb_records[[standard_idx[1]]]

all_mgi_mapping <- as.data.frame(genes(mouse_ensdb, columns = c("gene_id", "gene_name"))) %>%
  dplyr::rename(ensembl_gene_id = gene_id, mgi_symbol = gene_name) %>%
  dplyr::select(mgi_symbol, ensembl_gene_id) %>%
  
  mutate(mgi_symbol = as.character(mgi_symbol), ensembl_gene_id = as.character(ensembl_gene_id))

gene_symbols <- all_mgi_mapping %>%
  # Match on version-stripped Ensembl IDs as a defensive measure
 
  mutate(ensembl_gene_id_clean = sub("\\..*$", "", ensembl_gene_id)) %>%
  dplyr::filter(ensembl_gene_id_clean %in% sub("\\..*$", "", mouse_counts_ensembl$Geneid)) %>%
  dplyr::filter(!duplicated(ensembl_gene_id_clean)) %>%
  mutate(mgi_symbol = ifelse(is.na(mgi_symbol) | mgi_symbol == "", ensembl_gene_id, mgi_symbol)) %>%
  dplyr::select(mgi_symbol, ensembl_gene_id = ensembl_gene_id_clean)  # rename back so the join key downstream still works

# Genome-wide valid symbol list, reused later by the CMS section
# instead of re-querying Ensembl a second time.
valid_mgi_symbols <- all_mgi_mapping %>%
  dplyr::filter(!is.na(mgi_symbol), mgi_symbol != "") %>%
  pull(mgi_symbol) %>%
  unique()

mouse_counts_ensembl <- mouse_counts_ensembl %>%
  mutate(Geneid_clean = sub("\\..*$", "", Geneid)) %>%  # match gene_symbols' unversioned key; original Geneid column is preserved
  left_join(gene_symbols, by = c("Geneid_clean" = "ensembl_gene_id")) %>%
  dplyr::select(-Geneid_clean) %>%
  dplyr::filter(!is.na(mgi_symbol)) %>%
  mutate(mgi_symbol = make.unique(as.character(mgi_symbol)))
rownames(mouse_counts_ensembl) <- NULL

# ---------------------------------------------------------------
# 4. Raw counts matrix + gene lengths
# ---------------------------------------------------------------
message("4. Building counts matrix and gene lengths...")

mouse_counts_matrix <- mouse_counts_ensembl %>%
  dplyr::select(mgi_symbol, all_of(wanted_samples)) %>%
  mutate(across(all_of(wanted_samples), as.numeric)) %>%
  column_to_rownames("mgi_symbol") %>%
  as.matrix()

mouse_gene_lengths <- mouse_counts_ensembl %>%
  dplyr::select(mgi_symbol, length = Length) %>%
  mutate(length = as.numeric(length)) %>%
  column_to_rownames("mgi_symbol") %>%
  as.matrix()

# FINAL OUTPUT 1: TPM-normalized matrix
message("5. Converting to TPM...")
mouse_tpm_matrix <- convertCounts(mouse_counts_matrix, "TPM", geneLength = mouse_gene_lengths) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(!is.na(gene))
write.csv(mouse_tpm_matrix, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/TPMCounts.csv", row.names = TRUE)

# ---------------------------------------------------------------
# 6. Build dds_mouse
# ---------------------------------------------------------------
message("6. Building DESeq2 object from raw counts...")

mouse_counts_for_dds <- mouse_counts_ensembl %>%
  dplyr::select(-Geneid, -Length) %>%
  remove_rownames() %>%
  column_to_rownames("mgi_symbol")

mouse_metadata <- read.csv(metadata_file) %>%
  dplyr::rename(sample = sample_id_short) %>%
  mutate(
    sample    = gsub("\\.", "-", sample),
    bioSample = gsub("_Rep_[0-9]+$", "", sample)  # collapse technical reps
  )

count_order <- colnames(mouse_counts_for_dds)

mouse_metadata_aligned <- mouse_metadata
rownames(mouse_metadata_aligned) <- mouse_metadata_aligned$sample  # sets rownames WITHOUT removing "sample" -- collapseReplicates() needs it as a real column
mouse_metadata_aligned <- mouse_metadata_aligned[count_order, , drop = FALSE]

if(identical(rownames(mouse_metadata_aligned), colnames(mouse_counts_for_dds))) {
  print("Order Correct So will Proceed.")
}

dds_mouse <- DESeqDataSetFromMatrix(
  countData = mouse_counts_for_dds,
  colData   = mouse_metadata_aligned,
  design    = ~ grp
)
dds_mouse <- collapseReplicates(dds_mouse, dds_mouse$bioSample, dds_mouse$sample)
dds_mouse <- DESeq(dds_mouse)


write.csv(as.data.frame(assay(dds_mouse)),"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/GEO Submission/dds_all.csv", row.names=TRUE)

vsd_mouse <- vst(dds_mouse)
mouse_vst_matrix <- as.data.frame(assay(vsd_mouse))
mouse_vst_metadata <- as.data.frame(colData(dds_mouse)) %>%
  dplyr::select(-any_of("sample")) %>%   # stale technical-rep sample name, no longer accurate post-collapse
  rownames_to_column("sample")

write.csv(mouse_vst_matrix, file.path(data_output_dir, "mouse_vst_matrix.csv"), row.names = TRUE)
write.csv(mouse_vst_metadata, file.path(data_output_dir, "mouse_vst_metadata.csv"), row.names = FALSE)
message("Done building dds_mouse. Final objects: mouse_tpm_matrix, mouse_vst_matrix, mouse_vst_metadata.")

########################################################
# SUPPLEMENTARY UMAP -- all samples, colored by genotype,
# shaped by immune-exposure state
########################################################
umap_result <- umap(t(mouse_vst_matrix))
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df <- umap_df %>% rownames_to_column("sample") %>%
  merge(mouse_vst_metadata, by = "sample")

umap_df <- umap_df %>%
  mutate(
    genotype = case_when(
      grp %in% c("AK", "AK_Tu") ~ "AK",
      grp %in% c("AKP", "AKP_Tu", "AKP_TuOrganoid") ~ "AKP",
      grp %in% c("AKPS", "AKPS_Tu", "AKPS_TuOrganoid") ~ "AKPS_primary",
      grp == "NormalColon" ~ "NormalColon",
      grp == "Polyps" ~ "Polyps",
      grp == "AKPS_LungMet" ~ "AKPS_LungMet",
      grp == "AKPS_Met" ~ "AKPS_DistMet",
      grp == "AKPS_TdLN" ~ "AKPS_TdLN",
      TRUE ~ NA_character_  # flags any future/unexpected grp level instead of silently mis-plotting it
    ),
    immune_state = case_when(
      grp %in% c("AK", "AKP", "AKPS") ~ "Immune-naive organoid",
      grp %in% c("AKP_TuOrganoid", "AKPS_TuOrganoid") ~ "Immune-exposed organoid",
      grp %in% c("NormalColon", "Polyps") ~ "Normal/Polyp tissue",
      TRUE ~ "Spontaneous tumor"
    )
  )

genotype_colors <- c(
  AK = "#ff7f00", AKP = "#e41a1c", AKPS_primary = "brown",
  NormalColon = "gray40", Polyps = "#4daf4a",
  AKPS_LungMet = "violet", AKPS_DistMet = "#377eb8", AKPS_TdLN = "purple"
)
shape_values <- c(
  "Immune-naive organoid" = 16, "Immune-exposed organoid" = 17,
  "Normal/Polyp tissue" = 18, "Spontaneous tumor" = 15
)

write.csv(umap_df, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S2/umap.csv")

umap_Mouse_supp <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = genotype, shape = immune_state)) +
  geom_point(size = 5, stroke = 1) +
  scale_color_manual(values = genotype_colors) +
  scale_shape_manual(values = shape_values) +
  theme_classic(base_size = 15) +
  labs(x = "UMAP Dim 1", y = "UMAP Dim 2", color = "Genotype", shape = "Sample type")

ggsave(file.path(supp_fig_dir, "updatedUmap_Mouse.pdf"),
       plot = umap_Mouse_supp, width = 15, height = 6, units = "in", dpi = 300)

########################################################
# SUPPLEMENTARY correlation heatmap -- all samples 
########################################################
cor_matrix_all <- cor(mouse_vst_matrix, method = "spearman")
color_palette <- colorRampPalette(c("yellow", "orange", "red"))(100)

png(file.path(supp_fig_dir, "updatedHeatMapMouse.png"), width = 4000, height = 4000, res = 180)
corrplot(cor_matrix_all, tl.cex = 2, is.corr = FALSE, col = color_palette,
         cex.main = 3, cl.cex = 0.5, method = "color", type = "full",
         order = "hclust", hclust.method = "complete", title = "", mar = c(0, 0, 4, 0))
dev.off()

########################################################
# Volcano plot helpers (shared with main figures script)
########################################################
get_top_genes <- function(df, fc_filter, n = 10) {
  subset_df <- df[fc_filter(df$log2FoldChange), ]
  subset_df <- subset_df[order(subset_df$fdr_adj_pvalue), ]
  rownames(subset_df)[1:min(n, nrow(subset_df))]
}

make_volcano <- function(dds, contrast, plot_title, extra_labels = NULL,
                          fc_cut = 2, p_cut = 0.05, n_label = 10, show_legend = TRUE, printDataframe) {
  res <- lfcShrink(dds, contrast = contrast, type = "ashr")
  res$fdr_adj_pvalue <- p.adjust(res$pvalue, method = "fdr")

  up_genes   <- get_top_genes(res, function(fc) fc >= fc_cut, n = n_label)
  down_genes <- get_top_genes(res, function(fc) fc <= -fc_cut, n = n_label)
  select_labels <- unique(c(up_genes, down_genes, extra_labels, "Reg4"))

  p <- EnhancedVolcano(
    res, lab = rownames(res), selectLab = select_labels,
    x = "log2FoldChange", y = "fdr_adj_pvalue", title = plot_title,
    pCutoff = p_cut, FCcutoff = fc_cut, drawConnectors = TRUE,
    maxoverlapsConnectors = Inf, lengthConnectors = unit(0.0075, "npc"),
    labSize=20, 
    axisLabSize = 45,
    titleLabSize = 45,
    legendLabSize = 36,
    boxedLabels = FALSE
  )
  if (!show_legend) p <- p + theme(legend.position = "none")
  
  if(!is.null(printDataframe)) {
    write.xlsx(as.data.frame(res), paste0(printDataframe, "/", plot_title,".xlsx"), rowNames=TRUE)
  }
  
  p
}

########################################################
# SUPPLEMENTARY volcano: sequential progression (AK_Tu -> AKPS_Tu
# -> mets), all pulled from the full dds_mouse -- deliberate choice
# since this panel set specifically includes AK_Tu/mets, which are
# excluded from every main-figure dds object.
########################################################
volcano_seq_AKPTuvAKTu    <- make_volcano(dds_mouse, c("grp", "AKP_Tu", "AK_Tu"),        "AKP_Tu vs AK_Tu", show_legend = TRUE, printDataframe = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3")
volcano_seq_AKPSTuvAKPTu  <- make_volcano(dds_mouse, c("grp", "AKPS_Tu", "AKP_Tu"),      "AKPS_Tu vs. AKP_Tu", show_legend = FALSE, printDataframe = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3")
volcano_seq_LungMetvAKPSTu<- make_volcano(dds_mouse, c("grp", "AKPS_LungMet", "AKPS_Tu"),"LungMet vs AKPS_Tu", show_legend = FALSE, printDataframe = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3")
volcano_seq_TdLNvAKPSTu   <- make_volcano(dds_mouse, c("grp", "AKPS_TdLN", "AKPS_Tu"),   "TdLN vs AKPS_Tu", show_legend = FALSE, printDataframe = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3")
volcano_seq_MetvAKPSTu    <- make_volcano(dds_mouse, c("grp", "AKPS_Met", "AKPS_Tu"),    "DistMet vs AKPS_Tu", show_legend = FALSE, printDataframe = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3")

seq_legend <- get_legend(volcano_seq_AKPTuvAKTu + theme(legend.position = "bottom"))
volcano_seq_AKPTuvAKTu_clean <- volcano_seq_AKPTuvAKTu + theme(legend.position = "none")

seq_grid  <- plot_grid(volcano_seq_AKPTuvAKTu_clean, volcano_seq_AKPSTuvAKPTu, volcano_seq_LungMetvAKPSTu,
                        volcano_seq_TdLNvAKPSTu, volcano_seq_MetvAKPSTu, ncol = 3, align = "hv")
seq_final <- plot_grid(seq_grid, seq_legend, ncol = 1, rel_heights = c(1, 0.08))

pdf(file.path(supp_fig_dir, "VolcanoPlotSequentialMouse.pdf"), width = 16, height = 12)
tryCatch(print(seq_final), finally = dev.off())

########################################################
# SUPPLEMENTARY volcano: everything vs. Polyps (reviewer request)
########################################################
volcano_vp_AKTu   <- make_volcano(dds_mouse, c("grp", "AK_Tu", "Polyps"),   "AK Tumor vs Polyps", show_legend = TRUE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKPTu  <- make_volcano(dds_mouse, c("grp", "AKP_Tu", "Polyps"),  "AKP Tumor vs Polyps", show_legend = TRUE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKPSTu <- make_volcano(dds_mouse, c("grp", "AKPS_Tu", "Polyps"), "AKPS Tumor vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_LungMet<- make_volcano(dds_mouse, c("grp", "AKPS_LungMet", "Polyps"), "LungMet Tumor vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_TdLN   <- make_volcano(dds_mouse, c("grp", "AKPS_TdLN", "Polyps"),    "TdLN vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_Met    <- make_volcano(dds_mouse, c("grp", "AKPS_Met", "Polyps"),     "DistMet Tumor vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")


volcano_vp_AK     <- make_volcano(dds_mouse, c("grp", "AK", "Polyps"),   "AK organoid (immune naive) vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKP    <- make_volcano(dds_mouse, c("grp", "AKP", "Polyps"),  "AKP organoid (immune naive) vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKPS   <- make_volcano(dds_mouse, c("grp", "AKPS", "Polyps"), "AKPS organoid (immune naive) vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKPTuOrg  <- make_volcano(dds_mouse, c("grp", "AKP_TuOrganoid", "Polyps"),  "AKP organoid (immune exposed) vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")
volcano_vp_AKPSTuOrg <- make_volcano(dds_mouse, c("grp", "AKPS_TuOrganoid", "Polyps"), "AKPS organoid (immune exposed) vs Polyps", show_legend = FALSE, printDataframe="/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6")

vp_legend <- get_legend(volcano_vp_AKTu + theme(legend.position = "bottom"))
volcano_vp_AKTu_clean <- volcano_vp_AKTu + theme(legend.position = "none")

vp_grid  <- plot_grid(volcano_vp_AKTu_clean, volcano_vp_AKPTu, volcano_vp_AKPSTu, volcano_vp_LungMet,
                       volcano_vp_TdLN, volcano_vp_Met, volcano_vp_AK, volcano_vp_AKP, volcano_vp_AKPS,
                       volcano_vp_AKPTuOrg, volcano_vp_AKPSTuOrg, ncol = 5, align = "hv")
vp_final <- plot_grid(vp_grid, vp_legend, ncol = 1, rel_heights = c(1, 0.08))

pdf(file.path(supp_fig_dir, "VolcanoPlotAllvPolyp.pdf"), width = 35, height = 20)
tryCatch(print(vp_final), finally = dev.off())

write.csv(as.data.frame(assay(dds_mouse)), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S3/dds.csv")

########################################################
# SUPPLEMENTARY UpSet plot (Reviewer 3 suggestion)
########################################################
get_upregulated_genes <- function(dds, contrast, padj_cutoff = 0.05, lfc_cutoff = 0) {
  results(dds, contrast = contrast) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    dplyr::filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)
}

tumor_gene_lists_vs_polyps <- list(
  AK_Tumor      = get_upregulated_genes(dds_mouse, c("grp", "AK_Tu", "Polyps")),
  AKP_Tumor     = get_upregulated_genes(dds_mouse, c("grp", "AKP_Tu", "Polyps")),
  AKPS_Tumor    = get_upregulated_genes(dds_mouse, c("grp", "AKPS_Tu", "Polyps")),
  AK_Organoid   = get_upregulated_genes(dds_mouse, c("grp", "AK", "Polyps")),
  AKP_Organoid  = get_upregulated_genes(dds_mouse, c("grp", "AKP", "Polyps")),
  AKPS_Organoid = get_upregulated_genes(dds_mouse, c("grp", "AKPS", "Polyps"))
)

tumor_upset_df <- fromList(tumor_gene_lists_vs_polyps)

write.xlsx(tumor_upset_df,"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S6/TumorUpset.xlsx", rowNames=TRUE)

upset_plot <- ComplexUpset::upset(
  tumor_upset_df, names(tumor_gene_lists_vs_polyps),
  width_ratio = 0.2,
  base_annotations = list("Intersection size" = intersection_size(text = list(size = 3)))
) + theme(text = element_text(size = 12))

png(file.path(supp_fig_dir, "upsetPlotMouseVsPolyp.png"), width = 8000, height = 2000, res = 300)
print(upset_plot)
dev.off()

########################################################
# SUPPLEMENTARY immune cell stacked bar (8 disease stages)
########################################################
average_columns <- function(df, cols, new_name) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop("average_columns(): these columns are missing from the data: ",
         paste(missing_cols, collapse = ", "))
  }
  df[[new_name]] <- rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
  df[, !names(df) %in% setdiff(cols, new_name)]
}

collapse_technical_replicates <- function(df, id_col = "gene") {
  rep_pattern <- "(?i)_rep_[0-9]+$"
  sample_cols <- setdiff(names(df), id_col)
  tech_cols   <- grep(rep_pattern, sample_cols, value = TRUE, perl = TRUE)
  if (length(tech_cols) == 0) {
    message("No technical-replicate ('_Rep_N') columns found.")
    return(df)
  }
  base_names <- sub(rep_pattern, "", tech_cols, perl = TRUE)
  for (base in unique(base_names)) {
    group_cols <- tech_cols[base_names == base]
    df <- average_columns(df, group_cols, base)
  }
  df
}

mouse_tpm_averaged <- collapse_technical_replicates(mouse_tpm_matrix, id_col = "gene")


stage_definitions <- list(
  "Normal"       = c("NormalColon_002", "NormalColon_003", "NormalColon_004"),
  "Polyp"        = c("Polyp_1", "Polyps_2", "Polyps_3"),
  "AK_Tu"        = c("AK_Tu_1", "AK_Tu_2", "AK_Tu_3"),
  "AKP_Tu"       = c("AKP_Tu-1", "AKP_Tu-2"),
  "AKPS_Tu"      = c("AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3"),
  "AKPS_Met"     = c("AKPS_Met-1"),
  "AKPS_TdLN"    = c("AKPS_TdLN-1"),
  "AKPS_LungMet" = c("AKPS_LungMet-1")
)
for (stage in names(stage_definitions)) {
  mouse_tpm_averaged <- average_columns(mouse_tpm_averaged, stage_definitions[[stage]], stage)
}
mouse_tpm_averaged <- mouse_tpm_averaged %>% remove_rownames() %>% column_to_rownames("gene")

message("Running immune deconvolution (stage-averaged, supplementary)...")
res_deconv <- deconvolute_mouse(mouse_tpm_averaged, method = "mmcp_counter")

stage_order <- c("Normal", "Polyp", "AK_Tu", "AKP_Tu", "AKPS_Tu", "AKPS_Met", "AKPS_TdLN", "AKPS_LungMet")
res_deconv_matrix <- res_deconv %>% column_to_rownames("cell_type") %>% as.matrix()
res_deconv_matrix <- res_deconv_matrix[, stage_order]

cell_fractions <- res_deconv_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Cell_Type") %>%
  pivot_longer(cols = -Cell_Type, names_to = "Stage", values_to = "Fraction") %>%
  mutate(Cell_Type = case_when(
    Cell_Type %in% c("B cell", "B cell memory") ~ "B cell",
    Cell_Type %in% c("Monocyte", "Macrophage/Monocyte", "Granulocyte-monocyte progenitor") ~ "Myeloid Cells",
    TRUE ~ Cell_Type
  )) %>%
  dplyr::filter(!Cell_Type %in% c("Cancer associated fibroblast", "Endothelial cell")) %>%
  group_by(Stage, Cell_Type) %>%
  summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop") %>%
  mutate(Stage = factor(Stage, levels = stage_order))

stacked_bar_plot <- ggplot(cell_fractions, aes(x = Stage, y = Fraction, fill = Cell_Type)) +
  geom_col(position = "stack", alpha = 0.85, color = "white", linewidth = 0.2) +
  theme_minimal() +
  labs(title = "Cell Type Fractions Across Stages", y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "right")

write.xlsx(cell_fractions, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S5/cellFractions.xlsx")

pdf(file.path(supp_fig_dir, "StackedBarPlotMouse.pdf"))
tryCatch(print(stacked_bar_plot), finally = dev.off())

########################################################
# SUPPLEMENTARY GSEA: AK_Tu / AKP_Tu / AKPS_Tu / mets vs. Polyps
########################################################
baseline_grp <- "Polyps"
gsea_comparisons <- c(
  AK_Tu = "AK_Tu", AKP_Tu = "AKP_Tu", AKPS_Tu = "AKPS_Tu",
  AKPS_LungMet = "AKPS_LungMet", AKPS_TdLN = "AKPS_TdLN", AKPS_Met = "AKPS_Met"
)

make_preranked_file <- function(dds, grp_level, label, baseline, out_dir) {
  message(sprintf("  %s vs. %s...", label, baseline))
  res <- lfcShrink(dds, contrast = c("grp", grp_level, baseline), type = "ashr")
  ranked <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(ranking = sign(log2FoldChange) * -log10(pvalue)) %>%
    dplyr::select(gene, ranking) %>%
    dplyr::filter(is.finite(ranking))
  dup_genes <- ranked$gene[duplicated(ranked$gene)]
  if (length(dup_genes) > 0) {
    warning(sprintf("  %s: duplicate gene names in ranking output: %s",
                     label, paste(unique(dup_genes), collapse = ", ")))
  }
  ranked <- ranked[order(-ranked$ranking), ]
  out_path <- file.path(out_dir, sprintf("res_%svs%s.rnk", label, baseline))
  write.table(ranked, file = out_path, sep = "\t", row.names = FALSE,
              col.names = c("gene", ""), quote = FALSE)
  ranked
}

message("Building preranked .rnk files (vs. Polyps)...")
Map(make_preranked_file, grp_level = names(gsea_comparisons), label = unname(gsea_comparisons),
    MoreArgs = list(dds = dds_mouse, baseline = baseline_grp, out_dir = rnk_dir))

find_gsea_output_dir <- function(base_dir, label, baseline) {
  prefix  <- sprintf("%svs%s", label, baseline)
  matches <- list.dirs(base_dir, recursive = FALSE)
  matches <- matches[grepl(paste0("^", prefix, "\\.GseaPreranked\\."), basename(matches))]
  if (length(matches) == 0) {
    stop("No GSEA output folder found for comparison: ", prefix, " under ", base_dir,
         ". Has GSEA been run for this comparison yet?")
  }
  if (length(matches) > 1) {
    matches <- matches[order(file.info(matches)$mtime, decreasing = TRUE)]
    message(sprintf("    Multiple GSEA runs found for %s; using most recent.", prefix))
  }
  matches[1]
}

read_gsea_result <- function(label, baseline, base_dir, out_dir) {
  message(sprintf("  %s vs. %s...", label, baseline))
  run_dir  <- find_gsea_output_dir(base_dir, label, baseline)
  neg_file <- list.files(run_dir, pattern = "^gsea_report_for_na_neg_.*\\.tsv$", full.names = TRUE)
  pos_file <- list.files(run_dir, pattern = "^gsea_report_for_na_pos_.*\\.tsv$", full.names = TRUE)
  if (length(neg_file) != 1 || length(pos_file) != 1) {
    stop("Expected exactly one na_neg and one na_pos report in ", run_dir)
  }
  read_and_clean <- function(path) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(across(c(NES, `NOM p-val`, `FDR q-val`), as.numeric))
  }
  neg_report <- read_and_clean(neg_file)
  pos_report <- read_and_clean(pos_file)
  combined <- bind_rows(neg_report, pos_report) %>%
    dplyr::select(NAME, NES, `NOM p-val`, `FDR q-val`)
  write.xlsx(combined, file.path(out_dir, sprintf("res_%svs%s.xlsx", label, baseline)))
  combined
}

message("Reading and cleaning GSEA results (vs. Polyps)...")
gsea_results_supp <- Map(read_gsea_result, label = unname(gsea_comparisons),
                          MoreArgs = list(baseline = baseline_grp, base_dir = gsea_out_dir, out_dir = gsea_res_dir))
names(gsea_results_supp) <- unname(gsea_comparisons)

prepare_heatmap_data <- function(df, comparison_name) {
  df %>%
    dplyr::select(NAME, NES, `FDR q-val`) %>%
    mutate(Comparison = comparison_name, Significance = if_else(`FDR q-val` <= 0.05, "*", ""))
}

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

comparison_labels_supp <- sprintf("%s vs. %s", unname(gsea_comparisons), baseline_grp)
combined_data_supp <- Map(prepare_heatmap_data, gsea_results_supp, comparison_labels_supp) %>%
  bind_rows() %>%
  mutate(
    Pathway_Type = if_else(NAME %in% immune_pathways, "Immune", "Non-Immune"),
    Pathway_Order = if_else(Pathway_Type == "Immune", match(NAME, immune_pathways), match(NAME, oncogenic_pathways))
  ) %>%
  arrange(desc(Pathway_Type), Pathway_Order) %>%
  mutate(NAME = factor(NAME, levels = unique(NAME)))

combined_data_supp$Comparison <- factor(
  combined_data_supp$Comparison,
  levels = c("AK_Tu vs. Polyps", "AKP_Tu vs. Polyps", "AKPS_Tu vs. Polyps",
             "AKPS_Met vs. Polyps", "AKPS_TdLN vs. Polyps", "AKPS_LungMet vs. Polyps")
)

make_pathway_heatmap <- function(data, strip_text_size = 11) {
  ggplot(data, aes(Comparison, NAME, fill = NES)) +
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
}

ggsave(file.path(supp_fig_dir, "combined_pathway_heatmap_mouse_immune.pdf"),
       make_pathway_heatmap(combined_data_supp %>% filter(NAME %in% immune_pathways), 11), width = 18, height = 8)
ggsave(file.path(supp_fig_dir, "combined_pathway_heatmap_mouse_oncogenic.pdf"),
       make_pathway_heatmap(combined_data_supp %>% filter(NAME %in% oncogenic_pathways), 9), width = 18, height = 8)

########################################################
# SUPPLEMENTARY CMS classification (all 13 groups)
########################################################
message("Running MmCMS classification (supplementary, all groups)...")

mouse_vst_cms_input <- mouse_vst_matrix[rownames(mouse_vst_matrix) %in% valid_mgi_symbols, ]

cms_results <- MmCMS(
  as.matrix(mouse_vst_cms_input), templates = MmCMS::template.CMS.C,
  Genesets = c("template.CMS.C"), seed = 1, FDR = 1
)

cms_data <- cms_results %>%
  rownames_to_column("sample") %>%
  left_join(mouse_vst_metadata %>% dplyr::select(sample, grp), by = "sample") %>%
  mutate(
    prediction  = if_else(is.na(prediction), "Unclassified", prediction),
    group_label = gsub("_", " ", grp)
  )

group_order_all <- c("NormalColon", "Polyps", "AK_Tu", "AKP_Tu", "AKPS_Tu",
                      "AKPS_Met", "AKPS_LungMet", "AKPS_TdLN",
                      "AK", "AKP", "AKPS", "AKP_TuOrganoid", "AKPS_TuOrganoid")

cms_counts_all <- cms_data %>%
  dplyr::count(group_label = gsub("_", " ", grp), prediction) %>%
  mutate(
    group_label = factor(group_label, levels = gsub("_", " ", group_order_all)),
    prediction  = factor(prediction)
  )


write.xlsx(cms_data, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S5/CMSdata.xlsx")

stackedBar_Mouse_CMS_all <- ggplot(cms_counts_all, aes(x = group_label, y = n, fill = prediction)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1", na.value = "black") +
  labs(title = "CMS Subtype Distribution by Group (Counts)", x = "Group",
       y = "Number of Samples", fill = "CMS Subtype") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(supp_fig_dir, "stackedBar_Mouse_CMS.pdf"), plot = stackedBar_Mouse_CMS_all, width = 10, height = 8)

write.xlsx(mouse_vst_cms_input, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 6/mouse_vst_cms_input.xlsx")
write.xlsx(cms_data, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 6/cms_data.xlsx")


########## sup figure 8 ################### 
output_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Supplementary Figures"

# ---------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------
focal_groups <- c("Polyps", "AKP_Tu", "AKPS_Tu", "AK", "AKP", "AKPS")
top_n <- 200

# ---------------------------------------------------------------
# 1. Subset dds_mouse to just these 6 groups.
# ---------------------------------------------------------------
dds_six_groups <- dds_mouse[, dds_mouse$grp %in% focal_groups]
dds_six_groups$grp <- droplevels(dds_six_groups$grp)
message(sprintf("Subsetted to %d samples across %d groups: %s",
                ncol(dds_six_groups), length(unique(dds_six_groups$grp)),
                paste(levels(dds_six_groups$grp), collapse = ", ")))

run_one_vs_rest_top_genes <- function(dds, focal_group, top_n) {
  message(sprintf("  %s vs. rest...", focal_group))
  
  dds_binary <- dds
  colData(dds_binary)$one_vs_rest <- factor(
    ifelse(dds_binary$grp == focal_group, focal_group, "Other"),
    levels = c("Other", focal_group)
  )
  design(dds_binary) <- ~ one_vs_rest
  dds_binary <- DESeq(dds_binary)
  
  res <- lfcShrink(dds_binary, contrast = c("one_vs_rest", focal_group, "Other"), type = "ashr") %>%
    as.data.frame() %>%
    rownames_to_column("mouse_symbol") %>%
    dplyr::filter(padj < 0.05, log2FoldChange > 0) %>%   # "upregulated" = focal group > rest
    arrange(padj) %>%
    slice_head(n = top_n)
  
  message(sprintf("    %d genes passed the upregulated filter (requested top %d).", nrow(res), top_n))
  res$mouse_symbol
}

top_genes_by_group <- setNames(
  lapply(focal_groups, function(g) run_one_vs_rest_top_genes(dds_six_groups, g, top_n)),
  focal_groups
)

# ---------------------------------------------------------------
# 2. Map each group's top mouse genes to human Entrez orthologs
# ---------------------------------------------------------------
message("Mapping top genes to human orthologs...")

source("CrossSpeciesClassifier.R")

top_genes_human_entrez <- lapply(top_genes_by_group, function(genes) {
  mapped <- map_mouse_to_human_expr(genes, mouse_ensdb = mouse_ensdb)
  mapped <- mapped[!is.na(mapped$human_entrez), ]
  unique(mapped$human_entrez)
})

for (g in focal_groups) {
  message(sprintf("  %s: %d/%d top genes mapped to a human ortholog.",
                  g, length(top_genes_human_entrez[[g]]), length(top_genes_by_group[[g]])))
}

# ---------------------------------------------------------------
# 3. Overlap with CMScaller's CMS1-4 marker templates
# ---------------------------------------------------------------
cms_templates <- CMScaller::templates.CMS  # columns: probe (Entrez, character), class (CMS1-4)
cms_templates$probe <- as.character(cms_templates$probe)

overlap_df <- expand.grid(
  group = focal_groups,
  cms_class = levels(factor(cms_templates$class)),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_upregulated_mapped = length(top_genes_human_entrez[[group]]),
    n_overlap = length(intersect(
      top_genes_human_entrez[[group]],
      cms_templates$probe[cms_templates$class == cms_class]
    )),
    pct_overlap = if_else(n_upregulated_mapped > 0, 100 * n_overlap / n_upregulated_mapped, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    group = factor(group, levels = focal_groups),
    cms_class = factor(cms_class, levels = sort(unique(cms_class)))
  )

write.csv(overlap_df, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/S7/TopUnregulatedCMSOverlap.csv", row.names = FALSE)

# ---------------------------------------------------------------
# 4. Heatmap: x = group, y = CMS class, fill = % overlap (per PI suggestion)
# ---------------------------------------------------------------
overlap_heatmap <- ggplot(overlap_df, aes(x = group, y = cms_class, fill = pct_overlap)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f%%", pct_overlap)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "firebrick3", name = "% of top\nupregulated\ngenes", limits = c(0, NA)) +
  theme_minimal(base_size = 13) +
  labs(x = "Group (vs. rest)", y = "CMS class",
       title = sprintf("Overlap of top %d upregulated genes with CMS1-4 marker templates", top_n)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "top_upregulated_genes_CMS_overlap_heatmap.pdf"),
       overlap_heatmap, width = 9, height = 6)

message("Done.")
