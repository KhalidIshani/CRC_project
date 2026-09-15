

#Convert a counts matrix (genes x samples) + gene length to TPM
#counts_df data.frame with `width_kb` and one column per sample
#sample_cols character vector of sample column names

counts_to_tpm <- function(counts_df, sample_cols) {
  counts_df %>%
    mutate(across(all_of(sample_cols), ~ .x / width_kb)) %>%       # reads per kb
    mutate(across(all_of(sample_cols), ~ .x / sum(.x, na.rm = TRUE) * 1e6))  # scale to per-million
}



#---------------------------------------------------------------
# Shared helper: map row-name gene IDs to a new ID space, drop
# anything that fails to map, and collapse duplicates by mean.

# ---------------------------------------------------------------
# expr_mat matrix, genes (Ensembl IDs, versioned or not) x samples
#to_keytype target keytype for mapIds, e.g. "SYMBOL" or "ENTREZID"


map_and_collapse_ids <- function(expr_mat, to_keytype, label = to_keytype) {
  ensembl_clean <- sub("\\..*$", "", rownames(expr_mat))  # strip version suffix
  
  message(sprintf("Mapping Ensembl IDs to %s IDs...", label))
  mapped_ids <- mapIds(
    org.Hs.eg.db,
    keys      = ensembl_clean,
    column    = to_keytype,
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  
  is_mapped <- !is.na(mapped_ids)
  if (any(!is_mapped)) {
    warning(sprintf(
      "%d Ensembl IDs could not be mapped to %s IDs and were dropped.",
      sum(!is_mapped), label
    ))
  }
  
  expr_mapped <- expr_mat[is_mapped, ]
  ids_mapped  <- mapped_ids[is_mapped]
  
  if (any(duplicated(ids_mapped))) {
    message(sprintf("Averaging expression for duplicate %s IDs...", label))
    expr_mapped %>%
      as.data.frame() %>%
      mutate(.id = ids_mapped) %>%
      group_by(.id) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      column_to_rownames(".id")
  } else {
    out <- as.data.frame(expr_mapped)
    rownames(out) <- ids_mapped
    out
  }
}