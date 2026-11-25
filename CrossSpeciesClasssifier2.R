######################################################################
# Cross-species classifier (Figure 8) -- revised for R3C2
######################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "ensembldb", "DESeq2", "clusterProfiler", "org.Mm.eg.db"),
                     update = FALSE, ask = FALSE)
for (p in c("babelgene", "dplyr", "readr", "tibble", "nnet", "dynamicTreeCut", "msigdbr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(AnnotationHub)
library(ensembldb)
library(babelgene)
library(DESeq2)
library(dynamicTreeCut)
library(nnet)
library(dplyr)
library(readr)
library(tibble)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(clusterProfiler)

set.seed(100)

out_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis switches
R_THRESH            <- 0.6                     
R_THRESH_GRID       <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)
PROB_THRESH         <- 0.7                    
PRIMARY_TUMOR_ONLY  <- TRUE                    # keep TCGA "-01" samples only
RUN_STABILITY       <- TRUE                     # LOO + clustering null
N_CLUSTER_NULL      <- 10                      # permuted-data clusterings for the null

# ---------------------------------------------------------------
ah <- AnnotationHub()

mouse_ensdb_records <- query(ah, pattern = c("Mus musculus", "EnsDb", "110"))
print(mouse_ensdb_records)
mouse_ensdb <- mouse_ensdb_records[[1]]

# ---------------------------------------------------------------
# map_mouse_to_human_expr() 
# Returns a data.frame with one row per UNIQUE input mouse_symbol
map_mouse_to_human_expr <- function(mouse_symbols, mouse_ensdb, min_support = 3) {
  mouse_symbols <- unique(mouse_symbols)
  
  valid_genes <- as.data.frame(
    genes(mouse_ensdb, filter = ~ symbol %in% mouse_symbols, columns = "symbol")
  )
  valid_symbols <- unique(valid_genes$symbol)
  
  invalid <- setdiff(mouse_symbols, valid_symbols)
  if (length(invalid) > 0) {
    message(sprintf("  %d/%d input symbols not found in EnsDb v110 and will map to NA: %s",
                    length(invalid), length(mouse_symbols),
                    paste(head(invalid, 5), collapse = ", ")))
  }
  
  orth <- babelgene::orthologs(
    genes = valid_symbols, species = "mouse", human = FALSE,
    min_support = min_support, top = TRUE
  )
  
  mapped <- orth %>%
    dplyr::transmute(mouse_symbol = symbol, human_entrez = as.character(human_entrez)) %>%
    dplyr::distinct(mouse_symbol, .keep_all = TRUE)
  
  data.frame(mouse_symbol = mouse_symbols, stringsAsFactors = FALSE) %>%
    dplyr::left_join(mapped, by = "mouse_symbol")
}


######### cross-species classifier (immune-exposed tumors) ##########
tcga_zscore_file <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/coad_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"
comparing_samples <- c(
 "AK_1","AK_2","AK_3", "AKP_001", "AKP_004","AKP_005", "AKPS_1","AKPS_005","AKPS_006"
)

#comparing_samples <- c("Polyps_2","Polyps_3","Polyp_1","AKP_Tu-1","AKP_Tu-2_Rep_1", "AKP_Tu-2_Rep_2","AKPS_Tu-3","AKPS_Tu-2","AKPS_Tu-1")

mouse_counts_matrix <- read_rds("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/mouseCountsMatrix.rds")

cross_spec_counts <- mouse_counts_matrix[, comparing_samples]

sample_metadata <- data.frame(
  sample_id = colnames(cross_spec_counts),
  bio_id    = gsub("_Rep_[0-9]+$", "", colnames(cross_spec_counts)),
  group     = c(rep("AK", 3), rep("AKP", 3), rep("AKPS", 3))
)
rownames(sample_metadata) <- colnames(cross_spec_counts)

########################################################
# 2. DESeq2 -> VST normalization
########################################################
message("2. Collapsing replicates and computing VST...")

dds <- DESeqDataSetFromMatrix(
  countData = cross_spec_counts,
  colData   = sample_metadata,
  design    = ~ group
)

dds_collapsed <- collapseReplicates(dds, groupby = dds$bio_id, run = dds$sample_id)
dds_collapsed <- estimateSizeFactors(dds_collapsed)
vst_cross_spec <- vst(dds_collapsed, blind = TRUE)
cross_spec_norm <- assay(vst_cross_spec)  # genes x samples (bio_id-level)

bio_group_map <- unique(sample_metadata[, c("bio_id", "group")])

#######################################################
# 3. Top 8,000 most variable genes
########################################################
gene_vars <- apply(cross_spec_norm, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:8000]
expr_top  <- cross_spec_norm[top_genes, ]


########################################################
# 4. Hierarchical clustering (Spearman + Ward.D2)
########################################################
message("4. Clustering genes into modules...")

cluster_modules <- function(expr, min_size = 60) {
  cor_mat  <- cor(t(expr), method = "spearman")
  cor_mat[is.na(cor_mat)] <- 0
  dist_mat <- as.dist(1 - cor_mat)
  hc       <- hclust(dist_mat, method = "ward.D2")
  cl <- cutreeDynamic(
    dendro = hc, distM = as.matrix(dist_mat),
    minClusterSize = min_size, method = "hybrid", verbose = 0
  )
  names(cl) <- hc$labels
  gc <- split(names(cl), cl)
  gc[names(gc) != "0"]
}

gene_clusters <- cluster_modules(expr_top)

#######################################################
# 5. Module coherence filter 
########################################################
module_coherence <- function(gset, expr) {
  cs <- cor(t(expr[gset, , drop = FALSE]), method = "pearson")
  mean(cs[upper.tri(cs)], na.rm = TRUE)
}

# R3C2 (ii): fraction of variance explained by PC1 of the full mouse module
pc1_var_explained <- function(gset, expr) {
  x  <- scale(t(expr[gset, , drop = FALSE]))
  x  <- x[, colSums(is.na(x)) == 0, drop = FALSE]
  sv <- prcomp(x, center = FALSE, scale. = FALSE)$sdev^2
  sv[1] / sum(sv)
}

select_modules <- function(gene_clusters, expr, r_thresh, min_size = 60) {
  keep <- list()
  for (m in names(gene_clusters)) {
    gset <- gene_clusters[[m]]
    if (length(gset) < min_size) next
    if (module_coherence(gset, expr) > r_thresh) keep[[paste0("Module_", m)]] <- gset
  }
  keep
}

coherence_tbl <- data.frame(
  module         = paste0("Module_", names(gene_clusters)),
  n_genes_mouse  = lengths(gene_clusters),
  mean_within_r  = sapply(gene_clusters, module_coherence, expr = expr_top),
  mouse_pc1_var_full = sapply(gene_clusters, pc1_var_explained, expr = expr_top),
  stringsAsFactors = FALSE
)


write.csv(coherence_tbl, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model/CoherenceTblOrganoid.csv")



keep_modules <- select_modules(gene_clusters, expr_top, R_THRESH)
coherence_tbl$passes_r_thresh <- coherence_tbl$module %in% names(keep_modules)
message(sprintf("   %d modules retained with r > %.2f and >= 60 genes.", length(keep_modules), R_THRESH))


########################################################
# Shared gene map (computed once for all modules)
# R3C2: many-to-one orthologs collapsed so human rownames are unique
########################################################
message("Mapping mouse genes to human orthologs...")

gene_map_all <- map_mouse_to_human_expr(rownames(expr_top), mouse_ensdb = mouse_ensdb)
gene_map_all <- gene_map_all[!is.na(gene_map_all$human_entrez), ]

gene_map <- gene_map_all %>%
  dplyr::mutate(v = gene_vars[mouse_symbol]) %>%
  dplyr::arrange(dplyr::desc(v)) %>%
  dplyr::distinct(human_entrez, .keep_all = TRUE)

expr_top_human <- expr_top[gene_map$mouse_symbol, , drop = FALSE]
rownames(expr_top_human) <- gene_map$human_entrez

########################################################
# TCGA z-scores
########################################################
message("Reading TCGA z-scores...")

tcga_zscores <- read_delim(
  tcga_zscore_file, delim = "\t", show_col_types = FALSE,
  col_types = cols(Hugo_Symbol = col_character(),
                   Entrez_Gene_Id = col_character(),
                   .default = col_double())
) %>%
  dplyr::filter(!is.na(Entrez_Gene_Id), Entrez_Gene_Id != "") %>%
  dplyr::select(-dplyr::any_of("Hugo_Symbol")) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame() %>%
  column_to_rownames("Entrez_Gene_Id")

tcga_zscores <- as.matrix(tcga_zscores)

if (PRIMARY_TUMOR_ONLY) {
  tum <- grep("-01[A-Z]?$", colnames(tcga_zscores), value = TRUE)
  if (length(tum) > 0) {
    message(sprintf("   Keeping %d/%d primary-tumor samples.", length(tum), ncol(tcga_zscores)))
    tcga_zscores <- tcga_zscores[, tum, drop = FALSE]
  } else {
    message("   No '-01' sample IDs found; keeping all samples.")
  }
}

n_before <- nrow(tcga_zscores)
tcga_zscores <- tcga_zscores[complete.cases(tcga_zscores), ]
message(sprintf("   %d/%d genes retained after removing genes with missing values.",
                nrow(tcga_zscores), n_before))


########################################################
# R3C2 (A): shared-feature-space module scores
########################################################
project_module <- function(gset, mouse_expr, human_expr,
                           score_method = c("projected", "mean_z")) {
  score_method <- match.arg(score_method)
  gset <- intersect(gset, intersect(rownames(mouse_expr), rownames(human_expr)))
  if (length(gset) < 2) return(list(n_genes = length(gset)))
  
  m_s <- scale(t(mouse_expr[gset, , drop = FALSE]))   # samples x genes
  h_s <- scale(t(human_expr[gset, , drop = FALSE]))
  ok  <- colSums(is.na(m_s)) == 0 & colSums(is.na(h_s)) == 0
  m_s <- m_s[, ok, drop = FALSE]
  h_s <- h_s[, ok, drop = FALSE]
  if (ncol(m_s) < 2) return(list(n_genes = ncol(m_s)))
  
  pc_m <- prcomp(m_s, center = FALSE, scale. = FALSE)
  w    <- pc_m$rotation[, 1]
  if (cor(drop(m_s %*% w), rowMeans(m_s)) < 0) w <- -w
  
  pc_h  <- prcomp(h_s, center = FALSE, scale. = FALSE)   # diagnostic only
  h_axis <- drop(h_s %*% w)
  
  if (score_method == "projected") {
    mouse_score <- drop(m_s %*% w)
    human_score <- h_axis
  } else {
    mouse_score <- rowMeans(m_s)
    human_score <- rowMeans(h_s)
  }
  
  list(
    mouse   = mouse_score,
    human   = human_score,
    weights = w,
    n_genes = length(w),
    # R3C2 (ii)
    mouse_pc1_var_shared   = pc_m$sdev[1]^2 / sum(pc_m$sdev^2),
    human_var_by_mouse_axis = var(h_axis) / sum(apply(h_s, 2, var)),
    human_own_pc1_var      = pc_h$sdev[1]^2 / sum(pc_h$sdev^2),
    loading_cor_abs        = abs(cor(w, pc_h$rotation[, 1]))
  )
}

scale_safe <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

########################################################
# Steps 6-10 wrapped so they can be rerun for each threshold / scoring method
########################################################
run_classifier <- function(keep_modules, score_method = "projected",
                           min_genes = 10, prob_thresh = PROB_THRESH, verbose = TRUE) {
  stopifnot(all(grepl("^TCGA-", colnames(tcga_zscores))))
  # 6. Map mouse modules to human orthologs (>= 50% mapping, as before)
  human_gene_modules <- list()
  map_rate <- c()
  for (m in names(keep_modules)) {
    mouse_genes <- keep_modules[[m]]
    map_rate[m] <- mean(mouse_genes %in% gene_map_all$mouse_symbol)
    if (map_rate[m] >= 0.5) {
      human_gene_modules[[m]] <- gene_map$human_entrez[gene_map$mouse_symbol %in% mouse_genes]
    }
  }
  if (verbose) message(sprintf("   %d modules retained after >= 50%% mouse->human mapping.",
                               length(human_gene_modules)))
  
  # 7-8. Shared-space module scores (mouse = training, TCGA = prediction)
  proj <- lapply(human_gene_modules, project_module,
                 mouse_expr = expr_top_human, human_expr = tcga_zscores,
                 score_method = score_method)
  proj <- proj[vapply(proj, function(p) p$n_genes >= min_genes, logical(1))]
  if (length(proj) == 0) {
    warning("No modules left after mapping / shared-gene filtering.")
    return(NULL)
  }
  if (verbose) message(sprintf("   %d modules with >= %d genes shared with TCGA.", length(proj), min_genes))
  
  module_eigengenes <- as.data.frame(lapply(proj, `[[`, "mouse"))
  rownames(module_eigengenes) <- colnames(expr_top_human)
  
  eigengenes_human_mat <- as.data.frame(lapply(proj, `[[`, "human"))
  rownames(eigengenes_human_mat) <- colnames(tcga_zscores)
  
  diag_tbl <- data.frame(
    module                  = names(proj),
    mapping_rate            = unname(map_rate[names(proj)]),
    n_genes_shared          = sapply(proj, `[[`, "n_genes"),
    mouse_pc1_var_shared    = sapply(proj, `[[`, "mouse_pc1_var_shared"),
    human_var_by_mouse_axis = sapply(proj, `[[`, "human_var_by_mouse_axis"),
    human_own_pc1_var       = sapply(proj, `[[`, "human_own_pc1_var"),
    loading_cor_abs         = sapply(proj, `[[`, "loading_cor_abs"),
    stringsAsFactors = FALSE
  )
  
  # 9. Train multinomial model on mouse module scores 
  mouse_train <- module_eigengenes
  mouse_train$group <- bio_group_map$group[match(rownames(mouse_train), bio_group_map$bio_id)]
  stopifnot(!anyNA(mouse_train$group))
  
  predictor_cols <- setdiff(names(mouse_train), "group")
  mouse_train_scaled <- mouse_train
  mouse_train_scaled[predictor_cols] <- lapply(mouse_train_scaled[predictor_cols],
                                               function(x) as.numeric(scale(x)))
  
  multi_mod <- multinom(group ~ ., data = mouse_train_scaled, trace = FALSE, maxit = 1000)
  if (multi_mod$convergence != 0) warning("multinom did not converge (likely perfect separation).")
  
  # 10. Predict TCGA classes. All model features exist in TCGA by construction
  human_df <- eigengenes_human_mat[, predictor_cols, drop = FALSE]
  human_df_scaled <- as.data.frame(lapply(human_df, scale_safe))
  rownames(human_df_scaled) <- rownames(eigengenes_human_mat)
  
  human_pred_probs  <- predict(multi_mod, newdata = human_df_scaled, type = "probs")
  if (is.null(dim(human_pred_probs))) {                 # 1 TCGA sample edge case
    human_pred_probs <- matrix(human_pred_probs, nrow = 1,
                               dimnames = list(NULL, names(human_pred_probs)))
  }
  human_pred_labels <- as.character(predict(multi_mod, newdata = human_df_scaled, type = "class"))
  
  pred_results <- data.frame(
    Sample          = rownames(human_df_scaled),
    Predicted_Class = human_pred_labels,
    human_pred_probs,
    check.names = FALSE
  )
  prob_cols <- colnames(human_pred_probs)
  pred_results$highest_pred <- apply(pred_results[, prob_cols, drop = FALSE], 1, max)
  pred_results$Predicted_Class <- ifelse(
    pred_results$highest_pred < prob_thresh, "None", pred_results$Predicted_Class
  )
  
  list(pred_results = pred_results, diag = diag_tbl, model = multi_mod,
       mouse_scores = module_eigengenes, human_scores = eigengenes_human_mat)
}

########################################################
# Primary analysis
########################################################
message("6-10. Primary analysis (projected scores, r > ", R_THRESH, ")...")
main <- run_classifier(keep_modules, score_method = "projected")
print(summary(main$model))

write.csv(main$mouse_scores, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model", "module_scores_mouse_organoid.csv"))
write.csv(main$human_scores, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model", "module_scores_tcga_organoid.csv"))

##bind to clinical data 
main$pred_results$Sample <- substr(main$pred_results$Sample,1,15)

survival <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/Figure 8E-H/survivalData_organoid.csv")
names(survival)[names(survival) == "Predicted_Class"] <- "Predicted_Class_humanLoadings"

names(main$pred_results)[names(main$pred_results) == "Predicted_Class"] <- "Predicted_Class_mouseLoadings"

survivalData <- merge(main$pred_results, survival, by="Sample")

survival <- survival[survival$Sample %in% survivalData$Sample, ]


names(main$pred_results)[names(main$pred_results) == "Predicted_Class_mouseLoadings"] <- "Predicted_Class"






########################################################
# 11. R3C2 (i): GO / pathway enrichment per module
########################################################
message("11. Module enrichment (GO:BP, Hallmark)...")

out_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3"
  
go_list <- lapply(names(keep_modules), function(m) {
  eg <- tryCatch(enrichGO(
    gene = keep_modules[[m]], universe = top_genes,
    OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
    minGSSize = 10, maxGSSize = 500
  ), error = function(e) NULL)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
  eg <- clusterProfiler::simplify(eg, cutoff = 0.7, by = "p.adjust", select_fun = min)
  df <- as.data.frame(eg); df$module <- m; df
})
go_res <- dplyr::bind_rows(go_list)
write.csv(go_res, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model", "module_GO_BP_organoid_enrichment.csv"), row.names = FALSE)

# Hallmark on human orthologs (universe = mapped clustering genes)
hallmark <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                     error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
entrez_col <- intersect(c("ncbi_gene", "entrez_gene"), colnames(hallmark))[1]
hallmark_t2g <- data.frame(term = hallmark$gs_name, gene = as.character(hallmark[[entrez_col]]))

hm_list <- lapply(names(keep_modules), function(m) {
  hg <- gene_map$human_entrez[gene_map$mouse_symbol %in% keep_modules[[m]]]
  eh <- tryCatch(enricher(gene = hg, universe = gene_map$human_entrez,
                          TERM2GENE = hallmark_t2g, pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10),
                 error = function(e) NULL)
  if (is.null(eh) || nrow(as.data.frame(eh)) == 0) return(NULL)
  df <- as.data.frame(eh); df$module <- m; df
})
hm_res <- dplyr::bind_rows(hm_list)
write.csv(hm_res, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model", "module_Hallmark_enrichment_organoid.csv"), row.names = FALSE)

top_term <- function(res, m) {
  if (nrow(res) == 0 || !m %in% res$module) return(NA_character_)
  r <- res[res$module == m, ]; r$Description[which.min(r$p.adjust)]
}
n_terms <- function(res, m) if (nrow(res) == 0) 0L else sum(res$module == m)


########################################################
# 12. R3C2 (iii): justification + sensitivity of the r > 0.6 threshold
## (iii) Threshold justification + sensitivity
# 1. Null: coherence produced by the SAME clustering on structureless data
#    (each gene's values permuted independently across the 8 samples)
set.seed(1)
null_r <- unlist(replicate(50, {
  perm <- t(apply(expr_top, 1, sample)); dimnames(perm) <- dimnames(expr_top)
  gc <- cluster_modules(perm); gc <- gc[lengths(gc) >= 60]
  sapply(gc, module_coherence, expr = perm)
}, simplify = FALSE))
q99 <- quantile(null_r, 0.99)

coherence_tbl$passes_null99 <- coherence_tbl$mean_within_r > q99
print(q99)
print(coherence_tbl[, c("module", "mean_within_r", "passes_r_thresh", "passes_null99")])
identical(coherence_tbl$passes_r_thresh, coherence_tbl$passes_null99)   # want TRUE

# 2. Sensitivity: modules retained and agreement of TCGA calls with r > 0.6
grid <- t(sapply(c(0.5, 0.6, 0.7), function(th) {
  km  <- select_modules(gene_clusters, expr_top, th)
  res <- if (length(km)) run_classifier(km, verbose = FALSE) else NULL
  if (is.null(res)) return(c(r_thresh = th, n_modules = length(km), pct_agree = NA))
  res_ids <- substr(res$pred_results$Sample, 1, 15)
  calls <- res$pred_results$Predicted_Class[match(main$pred_results$Sample, res_ids)]
  c(r_thresh = th, n_modules = length(km), pct_agree = mean(calls == main$pred_results$Predicted_Class))
}))
print(grid)

write.csv(quantile(null_r, c(0.5, 0.95, 0.99)),"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model/NulldistOrganoid.csv")
write.csv(coherence_tbl[, c("module", "n_genes_mouse", "mean_within_r", "passes_r_thresh")],"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model/CoherenceTblOrganoid.csv")

write.csv(grid,"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/R3/Organoid Model/gridOrganoid.csv")
