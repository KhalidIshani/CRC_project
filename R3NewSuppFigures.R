######################################################################
# Supplementary Figures S11 (organoid) and S12 (immune-exposed tumor)
# Reviewer 3, Comment 2 -- built from saved files (no session rerun needed)
#
#   A  (iii) Module coherence vs. clustering null (95th / 99th pct) and r > 0.6
#   B  (ii)  Fraction of variance explained by PC1 per module
#   C  (i)   GO:BP enrichment (top non-redundant terms per retained module)
#   D  (i)   MSigDB Hallmark enrichment per retained module
#   E  R3C2  TCGA reclassification: original vs. corrected classifier
#
# PC1 variance: if the coherence CSV has no `mouse_pc1_var_full` column,
# the values in `pc1_var` below are used (taken from the canonical runs).
######################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(stringr); library(patchwork); library(readr)
})

######################################################################
# CONFIG -- set the two data folders to where the CSVs were written
######################################################################
base_dir     <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision"
organoid_dir <- file.path(base_dir, "DataSets/R3/Organoid Model")
tumor_dir    <- file.path(base_dir, "DataSets/R3/Full Tumor Model")
fig_dir      <- file.path(base_dir, "Supplementary Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
R_THRESH <- 0.6

models <- list(
  organoid = list(
    label        = "Supplementary Figure S11. Organoid-trained classifier: module validation",
    coherence    = file.path(organoid_dir, "CoherenceTblOrganoid.csv"),
    go           = file.path(organoid_dir, "module_GO_BP_organoid_enrichment.csv"),
    hallmark     = file.path(organoid_dir, "module_Hallmark_enrichment_organoid.csv"),
    null_q95     = 0.5624292,
    null_q99     = 0.5946110,
    pc1_var      = c(Module_1 = 0.6891021, Module_2 = 0.6300237, Module_3 = 0.6359859,
                     Module_4 = 0.7163944, Module_5 = 0.7163869, Module_6 = 0.7015078),
    reclass      = file.path(base_dir, "DataSets/Figure 8/Figure 8E-H/survivalData_organoid.csv"),
    new_col      = "Predicted_Class_mouseLoadings",
    old_col      = "Predicted_Class_humanLoadings_organoids",
    class_levels = c("AK", "AKP", "AKPS", "None"),
    out_file     = file.path(fig_dir, "R3SuppFiguresOrganoid.pdf")
  ),
  tumor = list(
    label        = "Supplementary Figure S12. Immune-exposed tumor classifier: module validation",
    coherence    = file.path(tumor_dir, "CoherenceTblTumor.csv"),
    go           = file.path(tumor_dir, "module_GO_BP_tumor_enrichment.csv"),
    hallmark     = file.path(tumor_dir, "module_Hallmark_enrichment_tumor.csv"),
    null_q95     = 0.6212441,
    null_q99     = 0.6455127,
    pc1_var      = c(Module_1 = 0.7582535, Module_2 = 0.7141048, Module_3 = 0.6519535,
                     Module_4 = 0.6640519, Module_5 = 0.5445035, Module_6 = 0.7089211),
    reclass      = file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8E-H/survivalData_TUMOR.csv"),  
    new_col      = "Predicted_Class_mouseLoadings",
    old_col      = "Predicted_Class_humanLoadings",
    class_levels = c("Polyps", "AKP_Tu", "AKPS_Tu", "None"),
    out_file     = file.path(fig_dir, "R3SuppFiguresTumor.pdf")
  )
)

######################################################################
# Helpers
######################################################################
status_cols <- c("Retained" = "#B2182B", "Excluded" = "#7F7F7F")
base_theme  <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

short_mod <- function(x) sub("^Module_", "M", x)

clean_hallmark <- function(x) {
  x <- str_to_sentence(gsub("_", " ", sub("^HALLMARK_", "", x)))
  for (a in c("E2f", "G2m", "Kras", "Myc", "Dna", "Uv", "Il2", "Il6", "Tnfa", "Jak",
              "Nfkb", "Stat3", "Stat5", "Pi3k", "Akt", "Mtorc1", "Tgf", "Dn", "Up", "V1", "V2"))
    x <- gsub(paste0("\\b", a, "\\b"), toupper(a), x, ignore.case = TRUE)
  x
}

gene_ratio <- function(gr) vapply(strsplit(gr, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 1)

# Top terms per module, skipping near-duplicates (gene-set Jaccard > max_jaccard)
top_nonredundant <- function(res, n_top = 4, max_jaccard = 0.5) {
  res %>% group_by(module) %>% arrange(p.adjust, .by_group = TRUE) %>%
    group_modify(function(d, ...) {
      keep <- integer(0); sets <- list()
      for (i in seq_len(nrow(d))) {
        g <- strsplit(d$geneID[i], "/")[[1]]
        redundant <- any(vapply(sets, function(s) length(intersect(s, g)) / length(union(s, g)) > max_jaccard, TRUE))
        if (!redundant) { keep <- c(keep, i); sets[[length(sets) + 1]] <- g }
        if (length(keep) == n_top) break
      }
      d[keep, ]
    }) %>% ungroup()
}

cohen_kappa <- function(a, b) {
  lv <- union(a, b); t <- table(factor(a, lv), factor(b, lv)); n <- sum(t)
  po <- sum(diag(t)) / n; pe <- sum(rowSums(t) * colSums(t)) / n^2
  (po - pe) / (1 - pe)
}

read_res <- function(path) {
  if (!file.exists(path)) { warning("Missing file: ", path); return(NULL) }
  r <- read_csv(path, show_col_types = FALSE)
  if (nrow(r) == 0) NULL else r
}

# Coherence table + PC1 variance (from CSV column if present, else from config)
read_coherence <- function(cfg, r_thresh) {
  ct <- as.data.frame(readr::read_csv(cfg$coherence, show_col_types = FALSE))
  message("Coherence file: ", cfg$coherence, "\n  columns: ", paste(names(ct), collapse = ", "))
  
  need <- c("module", "n_genes_mouse", "mean_within_r")
  if (!all(need %in% names(ct)))
    stop("Coherence CSV is missing: ", paste(setdiff(need, names(ct)), collapse = ", "))
  
  if (!"mouse_pc1_var_full" %in% names(ct)) {
    miss <- setdiff(ct$module, names(cfg$pc1_var))
    if (length(miss)) stop("No pc1_var in config for: ", paste(miss, collapse = ", "))
    ct$mouse_pc1_var_full <- unname(cfg$pc1_var[ct$module])
  }
  
  ct <- ct[, c(need, "mouse_pc1_var_full")]
  ct$status <- factor(ifelse(ct$mean_within_r > r_thresh, "Retained", "Excluded"),
                      levels = names(status_cols))
  ct$module <- factor(ct$module, levels = ct$module)
  ct
}

######################################################################
# Panels
######################################################################
panel_coherence <- function(ct, q95, q99, r_thresh) {
  ggplot(ct, aes(mean_within_r, module)) +
    annotate("rect", xmin = -Inf, xmax = q99, ymin = -Inf, ymax = Inf, fill = "grey93") +
    annotate("rect", xmin = -Inf, xmax = q95, ymin = -Inf, ymax = Inf, fill = "grey85") +
    geom_vline(xintercept = q95, linetype = "longdash", colour = "grey35") +
    geom_vline(xintercept = q99, linetype = "dashed") +
    geom_vline(xintercept = r_thresh, linetype = "dotted", colour = "#B2182B") +
    geom_segment(aes(x = min(c(0.45, mean_within_r)), xend = mean_within_r, yend = module),
                 colour = "grey70") +
    geom_point(aes(colour = status), size = 3) +
    geom_text(aes(label = sprintf("%.3f", mean_within_r)), hjust = -0.35, size = 2.9) +
    scale_colour_manual(values = status_cols, name = NULL) +
    scale_y_discrete(limits = rev(levels(ct$module)), labels = short_mod) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = "Module coherence vs. clustering null",
         subtitle = sprintf("Null 95th pct r = %.2f (long dash), 99th pct r = %.2f (dash)",
                            q95, q99),
         x = "Mean within-module Pearson r", y = NULL) +
    base_theme + theme(legend.position = "bottom")
}

panel_pc1 <- function(ct) {
  ggplot(ct, aes(module, mouse_pc1_var_full, fill = status)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * mouse_pc1_var_full)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = status_cols, name = NULL) +
    scale_x_discrete(labels = short_mod) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Variance explained by PC1",
         subtitle = "Full murine module",
         x = NULL, y = "Fraction of variance (PC1)") +
    base_theme + theme(legend.position = "bottom")
}

panel_enrich <- function(res, retained, title, clean = identity, n_top = 4, wrap = 38) {
  mods <- short_mod(retained)
  d <- if (is.null(res)) NULL else filter(res, module %in% retained)
  if (is.null(d) || nrow(d) == 0)
    return(ggplot() + annotate("text", 0, 0, label = "No significant terms") +
             theme_void() + labs(title = title))
  d <- d %>% top_nonredundant(n_top) %>%
    mutate(ratio = gene_ratio(GeneRatio),
           term  = str_wrap(clean(Description), wrap),
           mod   = factor(short_mod(module), levels = mods)) %>%
    arrange(mod, p.adjust)
  d$term <- factor(d$term, levels = rev(unique(d$term)))
  ns <- setdiff(mods, as.character(d$mod))
  ggplot(d, aes(mod, term)) +
    geom_point(aes(size = ratio, colour = -log10(p.adjust))) +
    scale_x_discrete(drop = FALSE) +
    scale_colour_gradient(low = "#F4A582", high = "#67001F", name = expression(-log[10]~adj.~p)) +
    scale_size_continuous(range = c(1.5, 6), name = "Gene ratio") +
    labs(title = title, x = "Module", y = NULL,
         subtitle = if (length(ns)) paste0("No significant terms: ", paste(ns, collapse = ", ")) else NULL) +
    base_theme + theme(axis.text.y = element_text(size = 7.5, lineheight = 0.85))
}

panel_reclass <- function(df, new_col, old_col, class_levels) {
  missing <- setdiff(c(new_col, old_col), names(df))
  if (length(missing))
    stop("Reclassification file lacks column(s): ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(names(df), collapse = ", "))
  new <- as.character(df[[new_col]]); old <- as.character(df[[old_col]])
  ok <- !is.na(new) & !is.na(old); new <- new[ok]; old <- old[ok]
  conf <- new != "None" & old != "None"
  tab <- as.data.frame(table(Corrected = factor(new, class_levels), Original = factor(old, class_levels)))
  ggplot(tab, aes(Original, Corrected, fill = Freq)) +
    geom_tile(colour = "white") +
    geom_tile(data = subset(tab, as.character(Original) == as.character(Corrected)),
              fill = NA, colour = "black", linewidth = 0.6) +
    geom_text(aes(label = Freq), size = 3.2) +
    scale_fill_gradient(low = "white", high = "#4393C3", guide = "none") +
    scale_y_discrete(limits = rev(class_levels)) +
    labs(title = "TCGA-COAD reclassification after correction",
         subtitle = sprintf("Agreement %d/%d (%.1f%%), kappa = %.2f; confident calls %.1f%%",
                            sum(new == old), length(new), 100 * mean(new == old),
                            cohen_kappa(new, old), 100 * mean(new[conf] == old[conf])),
         x = "Original (independent PC1)", y = "Corrected (fixed murine loadings)") +
    base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

######################################################################
# Figure
######################################################################
make_supp_module_figure <- function(cfg, r_thresh = R_THRESH, n_top = 4) {
  ct <- read_coherence(cfg, r_thresh)
  retained <- as.character(ct$module[ct$status == "Retained"])
  
  pA <- panel_coherence(ct, cfg$null_q95, cfg$null_q99, r_thresh)
  pB <- panel_pc1(ct)
  pC <- panel_enrich(read_res(cfg$go), retained, "GO Biological Process", n_top = n_top)
  pD <- panel_enrich(read_res(cfg$hallmark), retained, "MSigDB Hallmark",
                     clean = clean_hallmark, n_top = n_top)
  
  fig <- pA + pB + pC + pD +
    plot_layout(design = "AB\nCD", heights = c(1, 1.5)) +
    plot_annotation(title = cfg$label, tag_levels = "A",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(cfg$out_file, fig, width = 17, height = 12)
  message("Saved: ", cfg$out_file)
  invisible(fig)
}

for (m in models) make_supp_module_figure(m)