#===============================================================
# MOUSE COLON CANCER RNA-SEQ -- MAIN FIGURES
# Core group (NormalColon, Polyps, AKP_Tu, AKPS_Tu) and
# progenitor/organoid group (AKP, AKP_TuOrganoid, AKP_Tu, AKPS,
# AKPS_TuOrganoid, AKPS_Tu). Excludes AK_Tu and metastatic sites.
#
# REQUIRES 01_mouse_supplementary_figures.R to have been run
# first in the same session -- reuses dds_mouse, mouse_tpm_matrix,
# mouse_vst_matrix, mouse_vst_metadata, cms_data, valid_mgi_symbols,
# and helper functions
#===============================================================

main_fig_dir     <- file.path(project_dir, "Genome Medicine Revision/Main Figures")
main_fig_subdir  <- file.path(main_fig_dir, "Main Figures")
dir.create(main_fig_subdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# Subset + refit helper, reused throughout this script
# ---------------------------------------------------------------
subset_and_refit <- function(dds, groups) {
  dds_sub <- dds[, dds$grp %in% groups]
  dds_sub$grp <- droplevels(dds_sub$grp)
  message(sprintf("Subsetting to %d samples across %d groups: %s",
                   ncol(dds_sub), length(unique(dds_sub$grp)),
                   paste(levels(dds_sub$grp), collapse = ", ")))
  DESeq(dds_sub)
}

# ---------------------------------------------------------------
# Core group: NormalColon, Polyps, AKP_Tu, AKPS_Tu
# Progenitor/organoid group: AKP, AKPS, AKP_Tu, AKPS_Tu,
#   AKP_TuOrganoid, AKPS_TuOrganoid
# NOTE: AKPS_TuOrganoid is an n=1 biological-sample group after
# technical-replicate collapse
# ---------------------------------------------------------------
core_groups <- c("NormalColon", "Polyps", "AKP_Tu", "AKPS_Tu")
dds_core <- subset_and_refit(dds_mouse, core_groups)

progenitor_organoid_groups <- c("AKP", "AKPS", "AKP_Tu", "AKPS_Tu",
                                 "AKP_TuOrganoid", "AKPS_TuOrganoid")
dds_progenitor_organoid <- subset_and_refit(dds_mouse, progenitor_organoid_groups)

vst_core <- assay(vst(dds_core, blind = FALSE))
vst_progenitor_organoid <- assay(vst(dds_progenitor_organoid, blind = FALSE))

set.seed(100)
########################################################
# FIGURE 3C -- Core group UMAP
########################################################
umap_result_core <- umap(t(vst_core), n_neighbors = 6)  # reduced from default 15 -- too large for ~11 samples
umap_df_core <- as.data.frame(umap_result_core$layout)
colnames(umap_df_core) <- c("UMAP_1", "UMAP_2")
umap_df_core <- umap_df_core %>% rownames_to_column("sample") %>%
  merge(mouse_vst_metadata, by = "sample")

core_grp_colors <- c(AKP_Tu = "#e41a1c", NormalColon = "gray40", Polyps = "#4daf4a", AKPS_Tu = "#377eb8")

umap_Mouse_core <- ggplot(umap_df_core, aes(UMAP_1, UMAP_2, color = grp)) +
  geom_point(size = 5, stroke = 1) +
  scale_color_manual(values = core_grp_colors) +
  theme_classic(base_size = 15) +
  labs(x = "UMAP Dim 1", y = "UMAP Dim 2", color = "Group")

ggsave(file.path(main_fig_dir, "updatedUmap_Mouse_Core.pdf"),
       plot = umap_Mouse_core, width = 8, height = 6, units = "in", dpi = 300)

########################################################
# FIGURE 3A -- Progenitor/organoid group UMAP
########################################################
umap_result_po <- umap(t(vst_progenitor_organoid), n_neighbors = 9)  # ~14 samples
umap_df_po <- as.data.frame(umap_result_po$layout)
colnames(umap_df_po) <- c("UMAP_1", "UMAP_2")
umap_df_po <- umap_df_po %>% rownames_to_column("sample") %>%
  merge(mouse_vst_metadata, by = "sample")

po_grp_colors <- c(
  AKP_Tu = "#e41a1c", AKPS_Tu = "black", AKP = "#4daf4a",
  AKPS = "#377eb8", AKP_TuOrganoid = "orange", AKPS_TuOrganoid = "purple"
)

umap_Mouse_progenitor_organoid <- ggplot(umap_df_po, aes(UMAP_1, UMAP_2, color = grp)) +
  geom_point(size = 5, stroke = 1) +
  scale_color_manual(values = po_grp_colors) +
  theme_classic(base_size = 15) +
  labs(x = "UMAP Dim 1", y = "UMAP Dim 2", color = "Group")

ggsave(file.path(main_fig_dir, "updatedUmap_Mouse_Progenitor_Organoid.pdf"),
       plot = umap_Mouse_progenitor_organoid, width = 8, height = 6, units = "in", dpi = 300)

write.xlsx(umap_df_core, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3C/umap_df_core.xlsx", rowNames=TRUE)
write.xlsx(umap_df_po, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3A/umap_df_po.xlsx",rowNames=TRUE)
write.xlsx(vst_progenitor_organoid, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3A/vst_progenitor_organoid.xlsx",rowNames=TRUE)
write.xlsx(vst_core, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3C/vst_core.xlsx",rowNames=TRUE)

########################################################
# FIGURE 3B&D -- correlation heatmaps
########################################################
plot_corr_heatmap <- function(mat, out_path, clusters) {
  cor_matrix <- cor(mat, method = "spearman")
  color_palette <- colorRampPalette(c("yellow", "orange", "red"))(100)
  png(out_path, width = 4000, height = 4000, res = 180)
  corrplot(cor_matrix, tl.cex = 2, is.corr = FALSE, col = color_palette,
           cex.main = 3, cl.cex = 0.5, method = "color", type = "full",
           order = "hclust", hclust.method = "complete", addrect=clusters, title = "", mar = c(0, 0, 4, 0))
  dev.off()
}


# 3B immune-naive organoids only
immune_naive_organoid_groups <- c("AK", "AKP", "AKPS")
dds_immune_naive_organoid <- subset_and_refit(dds_mouse, immune_naive_organoid_groups)
vst_immune_naive_organoid <- assay(vst(dds_immune_naive_organoid, blind = FALSE))
plot_corr_heatmap(vst_immune_naive_organoid, file.path(main_fig_dir, "updatedHeatMapMouse_immune_naive_organoids.png"), 3)

# 3D core tumors only 
plot_corr_heatmap(vst_core, file.path(main_fig_dir, "updatedHeatMapMouse_core_tumors.png"), 4)


write.xlsx(vst_immune_naive_organoid, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3B/vst_immune_naive_organoid.xlsx",rowNames=TRUE)
write.xlsx(cor(vst_immune_naive_organoid, method="spearman"), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3B/correlations_immune_naive_organoid.xlsx",rowNames=TRUE)

write.xlsx(vst_core, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3D/vst_core.xlsx",rowNames=TRUE)
write.xlsx(cor(vst_core, method="spearman"), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3D/correlations_core.xlsx",rowNames=TRUE)


########################################################
# FIGURE 3E-L -- unified dds for all 8 volcano panels + GSEA
########################################################
groups_needed <- c("AK", "AKP", "AKPS", "AKP_Tu", "AKP_TuOrganoid",
                    "AKPS_Tu", "AKPS_TuOrganoid", "Polyps")
dds_all_included_samples <- subset_and_refit(dds_mouse, groups_needed)

path <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3E/"

volcano_AKPvAK <- make_volcano(
  dds_all_included_samples, c("grp", "AKP", "AK"),
  "AKP (immune-naive organoid) vs AK (immune-naive organoid)", show_legend = TRUE,
printDataframe=path)

volcano_AKPSvAKP <- make_volcano(
  dds_all_included_samples, c("grp", "AKPS", "AKP"),
  "AKPS (immune-naive organoid) vs AKP (immune-naive organoid)", show_legend = FALSE,printDataframe=path
)
volcano_AKPTuvAKP <- make_volcano(
  dds_all_included_samples, c("grp", "AKP_Tu", "AKP"),
  "AKP (spontaneous tumor) vs AKP (immune-naive organoid)", show_legend = FALSE,printDataframe=path
)
volcano_AKPTuOrganoidvAKP <- make_volcano(
  dds_all_included_samples, c("grp", "AKP_TuOrganoid", "AKP"),
  "AKP (immune-exposed tumor) vs AKP (immune-naive organoid)", show_legend = FALSE,printDataframe=path
)
volcano_AKPSTuvAKPS <- make_volcano(
  dds_all_included_samples, c("grp", "AKPS_Tu", "AKPS"),
  "AKPS (spontaneous tumor) vs AKPS (immune-naive organoid)", show_legend = FALSE,printDataframe=path
)
volcano_AKPSTuOrganoidvAKPS <- make_volcano(
  dds_all_included_samples, c("grp", "AKPS_TuOrganoid", "AKPS"),
  "AKPS (immune-exposed tumor) vs AKPS (immune-naive organoid)", show_legend = FALSE,printDataframe=path
)
volcano_AKPTuvPolyps <- make_volcano(
  dds_all_included_samples, c("grp", "AKP_Tu", "Polyps"),
  "AKP (spontaneous tumor) vs Polyps", show_legend = FALSE,printDataframe=path
)
volcano_AKPSTuvPolyps <- make_volcano(
  dds_all_included_samples, c("grp", "AKPS_Tu", "Polyps"),
  "AKPS (spontaneous tumor) vs Polyps", show_legend = FALSE,printDataframe=path
)

fig3c_legend <- get_legend(volcano_AKPvAK + theme(legend.position = "bottom"))
volcano_AKPvAK_clean <- volcano_AKPvAK + theme(legend.position = "none")

volcano_grid <- plot_grid(
  volcano_AKPvAK_clean, volcano_AKPTuvAKP, volcano_AKPSTuvAKPS, volcano_AKPTuvPolyps,
  volcano_AKPSvAKP, volcano_AKPTuOrganoidvAKP, volcano_AKPSTuOrganoidvAKPS, volcano_AKPSTuvPolyps,
  ncol = 4, align = "hv"
)
VolcanoPlotMainFigure <- plot_grid(volcano_grid, fig3c_legend, ncol = 1, rel_heights = c(1, 0.08))

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/Main Figures/Figure 3")
pdf("VolcanoPlotMainFigure.pdf", width = 95, height = 50)
VolcanoPlotMainFigure
dev.off()
write.xlsx(as.data.frame(assay(dds_all_included_samples)), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 3/Figure 3E/dds_all_included_samples.xlsx", rowNames=TRUE)


########################################################
# FIGURE 4 -- preranked GSEA
########################################################
make_preranked_file(dds_all_included_samples, "AKP_Tu",          "AKP_Tu",          "AKP",    rnk_dir)
make_preranked_file(dds_all_included_samples, "AKP_TuOrganoid",  "AKP_TuOrganoid",  "AKP",    rnk_dir)
make_preranked_file(dds_all_included_samples, "AKPS_Tu",         "AKPS_Tu",         "AKPS",   rnk_dir)
make_preranked_file(dds_all_included_samples, "AKPS_TuOrganoid", "AKPS_TuOrganoid", "AKPS",   rnk_dir)
make_preranked_file(dds_all_included_samples, "AKP_Tu",          "AKP_Tu",          "Polyps", rnk_dir)
make_preranked_file(dds_all_included_samples, "AKPS_Tu",         "AKPS_Tu",         "Polyps", rnk_dir)


make_preranked_file(dds_all_included_samples, "AKP_Tu",          "AKP_Tu",          "AKP_TuOrganoid", rnk_dir)
make_preranked_file(dds_all_included_samples, "AKPS_Tu",         "AKPS_Tu",         "AKPS_TuOrganoid", rnk_dir)




comparisons_df <- data.frame(
  label    = c("AKP_Tu", "AKP_TuOrganoid", "AKPS_Tu", "AKPS_TuOrganoid", "AKP_Tu", "AKPS_Tu", "AKP_Tu", "AKPS_Tu"),
  baseline = c("AKP",    "AKP",            "AKPS",    "AKPS",            "Polyps", "Polyps",   "AKP_TuOrganoid", "AKPS_TuOrganoid"),
  stringsAsFactors = FALSE
)

message("Reading and cleaning GSEA results (Figure 4A)...")
gsea_results <- Map(
  read_gsea_result,
  label    = comparisons_df$label,
  baseline = comparisons_df$baseline,
  MoreArgs = list(base_dir = gsea_out_dir, out_dir = gsea_res_dir)
)

names(gsea_results) <- comparisons_df$label
comparison_labels <- sprintf("%s vs. %s", comparisons_df$label, comparisons_df$baseline)

parent_idx <- c(1,2,3,4,7,8)   # AKP_Tu/AKP_TuOrganoid vs AKP, AKPS_Tu/AKPS_TuOrganoid vs AKPS
polyps_idx <- 5:6   # AKP_Tu vs Polyps, AKPS_Tu vs Polyps


build_pathway_heatmaps <- function(indices, file_prefix) {
  combined_data <- Map(prepare_heatmap_data, gsea_results[indices], comparison_labels[indices]) %>%
    bind_rows() %>%
    mutate(
      Pathway_Type = if_else(NAME %in% immune_pathways, "Immune", "Non-Immune"),
      Pathway_Order = if_else(Pathway_Type == "Immune", match(NAME, immune_pathways), match(NAME, oncogenic_pathways))
    ) %>%
    arrange(desc(Pathway_Type), Pathway_Order) %>%
    mutate(NAME = factor(NAME, levels = unique(NAME)))

  ggsave(file.path(main_fig_dir, sprintf("%s_immune.pdf", file_prefix)),
         make_pathway_heatmap(combined_data %>% filter(NAME %in% immune_pathways), 9),
         width = 18, height = 8)
  ggsave(file.path(main_fig_dir, sprintf("%s_oncogenic.pdf", file_prefix)),
         make_pathway_heatmap(combined_data %>% filter(NAME %in% oncogenic_pathways), 9),
         width = 18, height = 8)
  invisible(combined_data)
}

build_pathway_heatmaps(parent_idx, "tu_organoid_vs_parent_heatmap")
build_pathway_heatmaps(polyps_idx, "tu_vs_polyps_heatmap")

########################################################
# FIGURE 5 -- per-sample deconvolution + Wilcoxon comparisons
########################################################
message("Deconvolving per-sample TPM (Figure 5)...")

mouse_tpm_persample <- collapse_technical_replicates(mouse_tpm_matrix, id_col = "gene") %>%
  column_to_rownames("gene")
res_deconv_persample <- deconvolute_mouse(mouse_tpm_persample, method = "mmcp_counter")

deconv_long <- res_deconv_persample %>%
  pivot_longer(cols = -cell_type, names_to = "sample", values_to = "fraction") %>%
  left_join(mouse_vst_metadata %>% dplyr::select(sample, grp), by = "sample")

unmatched <- sum(is.na(deconv_long$grp))
if (unmatched > 0) warning(sprintf("%d rows had no matching grp in mouse_vst_metadata.", unmatched))

# NOTE: only B cell subtypes are merged here (not Monocyte/
# Macrophage/GM-progenitor into "Myeloid Cells" as in the
# supplementary stacked bar) 

cell_fractions_persample <- deconv_long %>%
  mutate(cell_type = case_when(
    cell_type %in% c("B cell", "B cell memory") ~ "B cell",
    TRUE ~ cell_type
  )) %>%
  dplyr::filter(!cell_type %in% c("Cancer associated fibroblast", "Endothelial cell")) %>%
  group_by(sample, grp, cell_type) %>%
  summarise(fraction = sum(fraction, na.rm = TRUE), .groups = "drop")

compare_two_groups <- function(long_data, group1, group2, alternative = "two.sided") {
  wide <- long_data %>%
    dplyr::filter(grp %in% c(group1, group2)) %>%
    dplyr::select(sample, grp, cell_type, fraction)
  cell_types <- unique(wide$cell_type)

  results <- lapply(cell_types, function(ct) {
    g1 <- wide$fraction[wide$cell_type == ct & wide$grp == group1]
    g2 <- wide$fraction[wide$cell_type == ct & wide$grp == group2]
    test <- tryCatch(
      wilcox.test(g1, g2, paired = FALSE, alternative = alternative),
      error = function(e) {
        warning(sprintf("Wilcoxon failed for %s (%s vs %s): %s", ct, group1, group2, e$message))
        list(p.value = NA)
      }
    )
    tibble(
      cell_type = ct,
      !!paste0("mean_", group1) := mean(g1, na.rm = TRUE),
      !!paste0("mean_", group2) := mean(g2, na.rm = TRUE),
      n_group1 = length(g1), n_group2 = length(g2),
      p_value = test$p.value
    )
  })

  bind_rows(results) %>%
    mutate(FDR_corrected_p_value = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

message("Running Figure 5 Wilcoxon comparisons...")
wilcox_comparisons <- list(c("NormalColon", "Polyps"), c("Polyps", "AKP"), c("AKP", "AKPS"))

for (pair in wilcox_comparisons) {
  res <- compare_two_groups(cell_fractions_persample, pair[1], pair[2])
  out_path <- file.path(main_fig_dir, sprintf("mouse_cellfractions_%svs%s.csv", pair[1], pair[2]))
  write.csv(res, out_path, row.names = FALSE)
  message(sprintf("   %s vs %s -> %s", pair[1], pair[2], out_path))
}

########################################################
# FIGURE 5 continued -- cell types split into separate
# panels/files by lineage, WITHOUT summing (each cell type keeps
# its own bar/color). Restricted to the core 4 groups.
########################################################
cell_panel_map <- c(
  "T cell" = "Lymphoid", "T cell CD8+" = "Lymphoid", "NK cell" = "Lymphoid",
  "B cell" = "Lymphoid", "B cell memory" = "Lymphoid",
  "Mast cell" = "Granulocytes", "Eosinophil" = "Granulocytes",
  "Neutrophil" = "Granulocytes", "Basophil" = "Granulocytes",
  "Monocyte" = "Myeloid", "Macrophage/Monocyte" = "Myeloid",
  "Granulocyte-monocyte progenitor" = "Myeloid"
)

plot_data <- cell_fractions_persample %>%
  dplyr::filter(cell_type %in% names(cell_panel_map)) %>%
  mutate(panel = cell_panel_map[cell_type]) %>%
  dplyr::filter(grp %in% core_groups) %>%
  mutate(grp = factor(grp, levels = core_groups)) %>%
  group_by(grp, panel, cell_type) %>%
  summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

for (pnl in unique(plot_data$panel)) {
  p <- ggplot(dplyr::filter(plot_data, panel == pnl), aes(x = grp, y = fraction, fill = cell_type)) +
    geom_col(position = "stack") +
    theme_minimal() +
    labs(x = "Group", y = "Mean Fraction", fill = "Cell Type", title = pnl) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(main_fig_dir, sprintf("immune_fractions_%s.pdf", gsub(" ", "_", pnl))),
         plot = p, width = 8, height = 6)
}

write.xlsx(plot_data, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 5/plotdata.xlsx")
write.xlsx(mouse_tpm_persample, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 5/mouse_tpm_persample.xlsx")

########################################################
# FIGURE 6 -- CMS 3-panel
########################################################
all_predictions <- c("CMS1", "CMS2", "CMS3", "CMS4")
stopifnot(all(setdiff(unique(cms_data$prediction), "Unclassified") %in% all_predictions))

pred_colors <- setNames(
  brewer.pal(max(3, length(all_predictions)), "Set1")[seq_along(all_predictions)],
  all_predictions
)

cms_counts <- cms_data %>%
  filter(prediction != "Unclassified") %>%
  dplyr::count(group_label = gsub("_", " ", grp), prediction) %>%
  tidyr::complete(group_label, prediction = all_predictions, fill = list(n = 0)) %>%
  mutate(prediction = factor(prediction, levels = all_predictions))

build_cms_barplot <- function(counts_df, group_order, title) {
  plot_data <- counts_df %>%
    filter(group_label %in% gsub("_", " ", group_order)) %>%
    mutate(group_label = factor(group_label, levels = gsub("_", " ", group_order)))

  ggplot(plot_data, aes(x = group_label, y = n, fill = prediction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pred_colors, drop = FALSE, na.value = "black") +
    labs(title = title, x = NULL, y = "Number of Samples", fill = "CMS Subtype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

akp_progression_order  <- c("AKP", "AKP_Tu", "AKP_TuOrganoid")
akps_progression_order <- c("AKPS", "AKPS_Tu", "AKPS_TuOrganoid")
full_group_order <- c("NormalColon", "Polyps", "AK_Tu", "AKP_Tu",
                       "AKPS_Tu", "AKPS_Met", "AKPS_LungMet", "AKPS_TdLN",
                       "AK", "AKP", "AKPS", "AKP_TuOrganoid", "AKPS_TuOrganoid")
core_group_order <- setdiff(full_group_order, c("AK_Tu", "AKPS_Met", "AKPS_LungMet", "AKPS_TdLN"))

cms_p1 <- build_cms_barplot(cms_counts, akp_progression_order,  "AKP Progression")
cms_p2 <- build_cms_barplot(cms_counts, akps_progression_order, "AKPS Progression")
cms_p3 <- build_cms_barplot(cms_counts, core_group_order,       "All Samples")


ggsave(file.path(main_fig_subdir, "CMS_AKP_progression.pdf"), plot = cms_p1)
ggsave(file.path(main_fig_subdir, "CMS_AKPS_progression.pdf"), plot = cms_p2)
ggsave(file.path(main_fig_subdir, "CMS_coreGroup_progression.pdf"), plot = cms_p3)

message("Main figures done.")
