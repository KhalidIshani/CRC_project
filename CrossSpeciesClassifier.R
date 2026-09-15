
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "ensembldb"), update = FALSE, ask = FALSE)
if (!requireNamespace("babelgene", quietly = TRUE)) install.packages("babelgene")

library(AnnotationHub)
library(ensembldb)
library(babelgene)


# ---------------------------------------------------------------
ah <- AnnotationHub()

mouse_ensdb_records <- query(ah, pattern = c("Mus musculus", "EnsDb", "110"))
print(mouse_ensdb_records)  
mouse_ensdb <- mouse_ensdb_records[[1]]

human_ensdb_records <- query(ah, pattern = c("Homo sapiens", "EnsDb", "110"))
print(human_ensdb_records)
human_mart <- human_ensdb_records[[1]] 

# ---------------------------------------------------------------
# map_mouse_to_human_expr()
#
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
  
  # Step 2: actual cross-species ortholog mapping 
  orth <- babelgene::orthologs(
    genes = valid_symbols, species = "mouse", human = FALSE,
    min_support = min_support, top = TRUE
  )
  
  mapped <- orth %>%
    dplyr::transmute(mouse_symbol = symbol, human_entrez = as.character(human_entrez)) %>%
    dplyr::distinct(mouse_symbol, .keep_all = TRUE)
  
  # Preserve every input gene as a row (NA where invalid or unmapped)
  data.frame(mouse_symbol = mouse_symbols, stringsAsFactors = FALSE) %>%
    dplyr::left_join(mapped, by = "mouse_symbol")
}


######### now build cross-species classifier ##########- this is going to be based on the immune-exposed organoids (not the spontaneous tumors)
tcga_zscore_file <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"
comparing_samples <- c(
  "Polyp_1", "Polyps_2", "Polyps_3",         
  "AKP_Tu-1", "AKP_Tu-2_Rep_1", "AKP_Tu-2_Rep_2",   
  "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3" 
)

cross_spec_counts <- mouse_counts_matrix[, comparing_samples]

#note that we will keep this general structure of separating sample_id and bio_id even though every sample in THIS CASE are biological replicates
sample_metadata <- data.frame(
  sample_id = colnames(cross_spec_counts),
  bio_id    = gsub("_Rep_[0-9]+$", "", colnames(cross_spec_counts)),
  group     = c(rep("Polyps", 3), rep("AKP_Tu", 3), rep("AKPS_Tu", 3))
)
rownames(sample_metadata) <- colnames(cross_spec_counts)

########################################################
# 2. DESeq2 -> VST normalization (no DE testing needed here)
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

cor_mat  <- cor(t(expr_top), method = "spearman")
dist_mat <- as.dist(1 - cor_mat)
hc       <- hclust(dist_mat, method = "ward.D2")

clusters <- cutreeDynamic(
  dendro = hc, distM = as.matrix(dist_mat),
  minClusterSize = 60, method = "hybrid"
)
names(clusters) <- hc$labels

gene_clusters <- split(names(clusters), clusters)
gene_clusters <- gene_clusters[names(gene_clusters) != "0"]

#######################################################
# 5. Keep modules with within-cluster correlation > 0.6
########################################################
keep_modules <- list()
for (m in names(gene_clusters)) {
  gset <- gene_clusters[[m]]
  if (length(gset) < 60) next
  cor_sub <- cor(t(expr_top[gset, ]), method = "pearson")
  avg_cor <- mean(cor_sub[upper.tri(cor_sub)])
  if (avg_cor > 0.6) keep_modules[[paste0("Module_", m)]] <- gset
}
message(sprintf("   %d modules retained with r > 0.6 and >= 60 genes.", length(keep_modules)))

########################################################
message("6. Mapping mouse gene modules to human orthologs...")

human_gene_modules <- list()
for (m in names(keep_modules)) {
  mouse_genes <- keep_modules[[m]]
  mapped <- map_mouse_to_human_expr(mouse_genes, mouse_ensdb = mouse_ensdb)
  mapped <- mapped[!is.na(mapped$human_entrez), ]
  
  mapping_rate <- nrow(mapped) / length(mouse_genes)
  if (mapping_rate >= 0.5) {
    human_gene_modules[[m]] <- as.character(mapped$human_entrez)
  }
}
message(sprintf("   %d modules retained after >= 50%% mouse->human mapping.", length(human_gene_modules)))

#######################################################
# 7. Mouse module eigengenes (training features)
########################################################
message("7. Computing mouse module eigengenes...")

mapped_expr <- map_mouse_to_human_expr(rownames(expr_top), mouse_ensdb = mouse_ensdb)

expr_top_human <- expr_top[rownames(expr_top) %in% mapped_expr$mouse_symbol, ]
rownames(expr_top_human) <- mapped_expr$human_entrez[
  match(rownames(expr_top_human), mapped_expr$mouse_symbol)
]

expr_top_human <- expr_top_human[!is.na(rownames(expr_top_human)), ]

module_eigengenes <- list()
for (m in names(human_gene_modules)) {
  gset <- intersect(human_gene_modules[[m]], rownames(expr_top_human))
  if (length(gset) >= 10) {
    module_data <- expr_top_human[gset, , drop = FALSE]
    eig <- WGCNA::moduleEigengenes(t(module_data), colors = rep(1, nrow(module_data)))$eigengenes[, 1]
    module_eigengenes[[m]] <- eig
  }
}
module_eigengenes <- as.data.frame(module_eigengenes)
rownames(module_eigengenes) <- colnames(expr_top_human)

########################################################
# 8. Human (TCGA) module eigengenes (prediction features)
########################################################
message("8. Computing human/TCGA module eigengenes...")

tcga_zscores <- read_delim(tcga_zscore_file, delim = "\t", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame() %>%
  column_to_rownames("Entrez_Gene_Id")

tcga_zscores <- as.matrix(tcga_zscores)
tcga_zscores <- tcga_zscores[complete.cases(tcga_zscores), ]

eigengenes_human_mat <- sapply(human_gene_modules, function(gset) {
  gset <- intersect(gset, rownames(tcga_zscores))
  if (length(gset) < 10) return(rep(NA, ncol(tcga_zscores)))
  sub_expr <- tcga_zscores[gset, , drop = FALSE]
  WGCNA::moduleEigengenes(t(sub_expr), colors = rep(1, nrow(sub_expr)))$eigengenes[, 1]
})
eigengenes_human_mat <- as.data.frame(eigengenes_human_mat)
rownames(eigengenes_human_mat) <- colnames(tcga_zscores)

########################################################
# 9. Train multinomial model on mouse module eigengenes
########################################################
message("9. Training model on mouse eigengenes...")

mouse_train <- module_eigengenes
mouse_train$group <- bio_group_map$group[match(rownames(mouse_train), bio_group_map$bio_id)]
stopifnot(!anyNA(mouse_train$group))

predictor_cols <- setdiff(names(mouse_train), "group")
mouse_train_scaled <- mouse_train
mouse_train_scaled[, predictor_cols] <- scale(mouse_train_scaled[, predictor_cols])

multi_mod <- multinom(group ~ ., data = mouse_train_scaled)
print(summary(multi_mod))

########################################################
# 10. Predict human sample labels from the mouse-trained model
########################################################
message("10. Predicting human sample classes...")

model_mods   <- predictor_cols
missing_mods <- setdiff(model_mods, colnames(eigengenes_human_mat))
if (length(missing_mods) > 0) {
  message(sprintf("   %d/%d model modules missing from human data; zero-imputing: %s",
                  length(missing_mods), length(model_mods), paste(missing_mods, collapse = ", ")))
}

human_df <- eigengenes_human_mat[, intersect(model_mods, colnames(eigengenes_human_mat)), drop = FALSE]
for (m in missing_mods) human_df[[m]] <- 0
human_df <- human_df[, model_mods, drop = FALSE]  # match model's column order exactly

# Scale each column independently on the human data's own
# distribution-- guarding against zero-variance columns
scale_safe <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}
human_df_scaled <- as.data.frame(lapply(human_df, scale_safe))
rownames(human_df_scaled) <- rownames(eigengenes_human_mat)

human_pred_probs  <- predict(multi_mod, newdata = human_df_scaled, type = "probs")
human_pred_labels <- as.character(predict(multi_mod, newdata = human_df_scaled, type = "class"))

pred_results <- data.frame(
  Sample          = rownames(human_df_scaled),
  Predicted_Class = human_pred_labels,
  human_pred_probs,
  check.names = FALSE
)

prob_cols <- colnames(human_pred_probs)  # actual class-probability columns, not hardcoded positions
pred_results$highest_pred <- apply(pred_results[, prob_cols, drop = FALSE], 1, max)


pred_results$Predicted_Class <- ifelse(
  pred_results$highest_pred < 0.7, "None", pred_results$Predicted_Class
)

write.csv(pred_results, file.path("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision/DataSets/Figure 8/", "all_pred_results_TumorResults.csv"), row.names = FALSE)

message("Done.")
