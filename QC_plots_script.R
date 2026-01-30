library(tidyverse)
library(ggplot2)
library(umap)
library(UpSetR)
library(RColorBrewer)
library(corrplot)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(DGEobj.utils)
library(janitor)
library(biomaRt)
library(tidyverse)
library(ggalluvial)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(gridExtra)
install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")
library(immunedeconv)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
library(networkD3)
library(ggpubr)
library(cowplot)
library(TCGAbiolinks)
library(SummarizedExperiment) # Used for normalization in some contexts, but not directly in this TPM calc.
library(immunedeconv) # For CIBERSORT (or other deconvolution methods)
library(org.Hs.eg.db) # For gene ID mapping
library(CMSclassifier) # For human CMS classification
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(scales)
library(ggalluvial)


##########now make new dataframe for immunedeconv ############
exp54 <- as.data.frame(data.table::fread("/Users/khalidishani/Downloads/subread_counts-2.txt", sep="\t"))
exp51 <- as.data.frame(data.table::fread("/Users/khalidishani/Downloads/subread_counts.txt", sep="\t"))
totalexp <- bind_cols(exp54,exp51)

wanted_samples <- c("AK_Tu_1", "AK_Tu_2","AK_Tu_3","Normal_002_1", "Normal_002_2","Normal_002_3", "Polyps_2", "Polyps_3",            
"Polyp_1","AKPS_LungMet-1_Rep_1", "AKPS_LungMet-1_Rep_2", "AKPS_LungMet-1_Rep_3",
"AKPS_Met-1_Rep_1","AKPS_Met-1_Rep_2","AKPS_Met-1_Rep_3","AKPS_TdLN-1_Rep_1",   
"AKPS_TdLN-1_Rep_2","AKPS_TdLN-1_Rep_3","AKPS_Tu-1","AKPS_Tu-2",           
"AKPS_Tu-3","AKP_Tu-1","AKP_Tu-2_Rep_1","AKP_Tu-2_Rep_2")

dmmr_pole_counts_ensembl <- totalexp[,c(wanted_samples,"Geneid...1","Length...6")]

dmmr_pole_counts_ensembl <- dmmr_pole_counts_ensembl %>% rename(
    Geneid...1 = "Geneid", 
    Length...6 = "Length"
)

#this dmmr_pole_counts_ensembl df is the unnormalized counts 
dmmr_pole_counts_ensembl <- dmmr_pole_counts_ensembl[rowSums(dmmr_pole_counts_ensembl==0) < ncol(dmmr_pole_counts_ensembl) -2, ]

# Get mapping of ensembl gene ids to gene symbols
# If ensembl id matches to "" empty string, replace with ensembl id instead

ensembl_106 <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl',version = 110, host="https://www.ensembl.org")

gene_symbols <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'), mart = ensembl_106) %>%
  filter(ensembl_gene_id %in% dmmr_pole_counts_ensembl$Geneid) %>%
  filter(!duplicated(ensembl_gene_id)) %>%
  mutate(mgi_symbol = ifelse(mgi_symbol == "", ensembl_gene_id, mgi_symbol))

# Add gene lengths to ensembl version of counts and make gene names unique
dmmr_pole_counts_ensembl <- dmmr_pole_counts_ensembl %>%
 #dplyr::select(-mgi_symbol) %>%
  left_join(., gene_symbols, by = c("Geneid" = "ensembl_gene_id")) %>%
  mutate(mgi_symbol = make.unique(mgi_symbol))


dmmr_pole_counts_ensembl <- dmmr_pole_counts_ensembl[!is.na(dmmr_pole_counts_ensembl$mgi_symbol),]

rownames(dmmr_pole_counts_ensembl) <- NULL

# Extract matrix of just the counts
dmmr_pole_counts <- dmmr_pole_counts_ensembl %>%
  dplyr::select(mgi_symbol, all_of(wanted_samples)) %>%
  mutate(across(all_of(wanted_samples),as.numeric)) %>%
  column_to_rownames("mgi_symbol") %>%
  as.matrix()

# Extract matrix of just the gene lengths
dmmr_pole_lengths <- dmmr_pole_counts_ensembl %>%
  dplyr::select(mgi_symbol, length = Length) %>%
  mutate(length = as.numeric(length)) %>%
  column_to_rownames("mgi_symbol") %>%
  as.matrix()

# Convert to TPM
#this is TPM normalized data frame 
dmmr_pole_tpm <- convertCounts(dmmr_pole_counts, "TPM", geneLength = dmmr_pole_lengths) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(gene))

################NOW DO QUALITY CONTROL USING RAW COUNTS -(df = dmmr_pole_counts_ensembl)- ###########################
setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output")

dmmr_pole_counts_ensembl_dds <- dmmr_pole_counts_ensembl[,!colnames(dmmr_pole_counts_ensembl) %in% c("Length","Geneid")]
dmmr_pole_counts_ensembl_dds<- dmmr_pole_counts_ensembl_dds %>% remove_rownames %>% column_to_rownames(var="mgi_symbol")
colonCancerMetaData <-  read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/all_samples_meta.csv")
colonCancerMetaData <- colonCancerMetaData[colonCancerMetaData$vitro_vivo == "in_vivo",]
colonCancerMetaData <- colonCancerMetaData[!colonCancerMetaData$grp %in% c("NormalColon","TdLN"),]
names(colonCancerMetaData)[names(colonCancerMetaData)=="sample_id_short"] <- "sample"
colonCancerMetaData$sample <- gsub("\\.", "-", colonCancerMetaData$sample)

#use VST method to normalize data 

##here would be fix 

colonCancerMetaData <- colonCancerMetaData %>%
  mutate(bioSample = gsub("_Rep_[0-9]+$", "", sample)) %>% 
  mutate(bioSample = make.unique(bioSample)) # remove _Rep_1, _Rep_2 from technical reps


dds_colonCancer <- DESeqDataSetFromMatrix(countData = dmmr_pole_counts_ensembl_dds,
                                          colData = colonCancerMetaData ,
                                          design= ~grp)


dds_colonCancer <- collapseReplicates(dds_colonCancer, dds_colonCancer$bioSample, dds_colonCancer$sample)


dds_colonCancer <- DESeq(dds_colonCancer)
vsd_colonCancer <- vst(dds_colonCancer)

data <- as.data.frame(assay(vsd_colonCancer))

colonCancerMetaData <- as.data.frame(colData(dds_colonCancer))

colonCancerMetaData <- colonCancerMetaData[,-1]

colonCancerMetaData <- colonCancerMetaData %>% rownames_to_column("sample")



setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output")
# =========================================================================
# Generate All 5 QC Plots

# -------------------------------------------------------------------------
# PCA

set.seed(123)

pca_result <- prcomp(t(data))
pca_df <- as.data.frame(pca_result$x)
pca_df <- pca_df %>%
  rownames_to_column("sample")
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
#KI addition - add coloring based on sample group 
pca_df <- merge(pca_df, colonCancerMetaData, by="sample")
#end KI addition
pdf("pca_plot.pdf", width = 8, height = 9)
ggplot(pca_df, aes(x = PC1, y = PC2, color=grp)) +
  geom_point(size=7) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(variance_explained[1], 2), "% variance)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "% variance)"))  +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=0),
        axis.ticks = element_line(linewidth = 0.5)) 
dev.off()

# -------------------------------------------------------------------------
# UMAP

umap_result <- umap(t(data))

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df <- umap_df %>%
  rownames_to_column("sample")

#KI addition - add coloring based on sample group 
umap_df <- merge(umap_df, colonCancerMetaData, by="sample")
#end KI addition

png("umap_plot_Nov6.png", width = 9, height = 9)
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color=grp)) +
  geom_point(size=6) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 0),
        axis.ticks = element_line(size = 0.5)) +
  labs(title = "UMAP Plot", x = "UMAP Dim 1", y = "UMAP Dim 2")
dev.off()


# -------------------------------------------------------------------------
# Multi-Line Density Plot

long_data <- data  %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

#KI addition - add coloring based on sample group 
long_data <- merge(long_data, colonCancerMetaData, by="sample")
#end KI addition

pdf("multi_line_density_plot.pdf", width = 10, height = 7)
ggplot(long_data, aes(x = expression, group = sample, color=grp)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Multi-Line Density Plot", x ="Expression", y = "Density") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),axis.line = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.ticks = element_line(size = 0.5)) 
dev.off()

# -------------------------------------------------------------------------
# Multi-Line Cutoff Plot

results_list <- c()
cutoffs <- seq(0,max(data ),0.25)

for (i in 1:length(cutoffs)) {
  c <- cutoffs[i]
  results_list[[i]] <- colSums(data > c)
}

num_genes_df <- as.data.frame(t(matrix(unlist(results_list), nrow=ncol(data),ncol=length(cutoffs))))
colnames(num_genes_df) <- colnames(data)
num_genes_df$Cutoff <- cutoffs
num_genes_df <- num_genes_df %>% 
  dplyr::select(Cutoff, everything()) %>%
  pivot_longer(-Cutoff, names_to = "sample", values_to = "NumberOfGenes") 

num_genes_colored <- num_genes_df %>%
  group_by(sample) %>%
  mutate(trend_y = mean(NumberOfGenes)) %>%
  ungroup() %>%
  mutate(total_deviation = abs(NumberOfGenes - trend_y)) %>%
  group_by(sample) %>%
  summarise(deviation = mean(total_deviation, na.rm = TRUE), 
            .groups = 'drop') %>%
  mutate(
    mean_dev = mean(deviation, na.rm = TRUE),
    sd_dev = sd(deviation, na.rm = TRUE),
    threshold = mean_dev + 2 * ifelse(is.na(sd_dev), 0, sd_dev),  # Handle NA sd
    potential_outlier = ifelse(deviation > threshold | deviation < mean_dev - 2.4 * ifelse(is.na(sd_dev), 0, sd_dev), "Outlier", "Normal"), # May need to play around with this to get the right visual outlier lines colored
    outlier_color = ifelse(potential_outlier == "Outlier", "pink", "black")
  ) %>%
  left_join(num_genes_df, by = "sample") %>%
  dplyr::select(Cutoff, sample, NumberOfGenes, outlier_color)

#KI addition - add coloring based on sample group 
num_genes_colored <- merge(num_genes_colored, colonCancerMetaData, by="sample")
#end KI addition

pdf("cutoff_plot.pdf",width = 10, height = 7)
ggplot(num_genes_colored, aes(x = Cutoff, y = NumberOfGenes, group = sample, color=grp)) +
  geom_line() +
  labs(title = "Cutoff Plot",
       x ="Expression",
       y = "Number of Genes > Cutoff",
       color = "sample") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),axis.line = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),
        axis.ticks = element_line(size = 0.5))
dev.off()


# -------------------------------------------------------------------------
# Sample-to-Sample ALAN2 UMAP Plot

cor_matrix_alan1 <- cor((data), method = "spearman")
cor_matrix_alan2 <- cor((cor_matrix_alan1), method = "pearson")
alan2_umap_result <- umap(cor_matrix_alan2)

alan2_umap_df <- as.data.frame(alan2_umap_result$layout)
colnames(alan2_umap_df) <- c("UMAP_1", "UMAP_2")
alan2_umap_df <- alan2_umap_df %>%
  rownames_to_column("sample")

#KI addition - add coloring based on sample group 
#alan2_umap_df$sample <- rownames(alan2_umap_df)
alan2_umap_df <- merge(alan2_umap_df, colonCancerMetaData, by="sample")
alan2_umap_df <- column_to_rownames(alan2_umap_df,"sample")
#end KI addition

pdf("alan2_umap_plot.pdf", width = 9, height = 9)
ggplot(alan2_umap_df, aes(x = UMAP_1, y = UMAP_2, color=grp)) +
  geom_point() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),
        axis.ticks = element_line(size = 0.5)) + 
  labs(title = "Sample-to-Sample ALAN2 UMAP", x = "UMAP 1", y = "UMAP 2")  + 
  theme(legend.position = "right") 
dev.off()


###make sample-to-sample covariance plot 
hclust_order <- hclust(dist(cor_matrix_alan1))
ordered_colnames <- colnames(cor_matrix_alan1)[order.dendrogram(as.dendrogram(hclust_order))]


color_palette <- colorRampPalette(c("yellow", "orange", "red"))(100)

png("cutoff_plot.png", width = 3000, height = 2000, res = 180)

# Create the heatmap
  corrplot(cor_matrix_alan1, 
           tl.cex=2,
         is.corr = FALSE,
         col = color_palette,  # Use the custom color palette
         cex.main = 3,            
         cl.cex = 0.5, 
         method = "color",          
         type = "full",            
         order = "hclust",         
         hclust.method = 'complete', 
         title = "",
         mar = c(0, 0, 4, 0))

dev.off()


########################################now do all averaging stuff - WITH TPM ##################################
dmmr_pole_tpm_AVERAGED <- dmmr_pole_tpm

dmmr_pole_tpm_AVERAGED$AKPS_LungMet <- (dmmr_pole_tpm_AVERAGED$`AKPS_LungMet-1_Rep_1` +
                                   dmmr_pole_tpm_AVERAGED$`AKPS_LungMet-1_Rep_2` +
                                   dmmr_pole_tpm_AVERAGED$`AKPS_LungMet-1_Rep_3`) / 3

dmmr_pole_tpm_AVERAGED$AKPS_Met <- (dmmr_pole_tpm_AVERAGED$`AKPS_Met-1_Rep_1` +
                               dmmr_pole_tpm_AVERAGED$`AKPS_Met-1_Rep_2` +
                               dmmr_pole_tpm_AVERAGED$`AKPS_Met-1_Rep_3`) / 3
dmmr_pole_tpm_AVERAGED$AKPS_TdLN <- (dmmr_pole_tpm_AVERAGED$`AKPS_TdLN-1_Rep_1` +
                                dmmr_pole_tpm_AVERAGED$`AKPS_TdLN-1_Rep_2` +
                                dmmr_pole_tpm_AVERAGED$`AKPS_TdLN-1_Rep_3`) / 3

dmmr_pole_tpm_AVERAGED$`AKP_Tu-2` <- (dmmr_pole_tpm_AVERAGED$`AKP_Tu-2_Rep_1` +
                                 dmmr_pole_tpm_AVERAGED$`AKP_Tu-2_Rep_2`) / 2

dmmr_pole_tpm_AVERAGED$Normal_002 <- (dmmr_pole_tpm_AVERAGED$Normal_002_1 + dmmr_pole_tpm_AVERAGED$Normal_002_2 + dmmr_pole_tpm_AVERAGED$Normal_002_3)/3


#this is the list of column names that were treated as technical replicates and averaged
dmmr_pole_tpm_AVERAGED <- dmmr_pole_tpm_AVERAGED %>%
  dplyr::select(-c(`AKPS_LungMet-1_Rep_1`, `AKPS_LungMet-1_Rep_2`, `AKPS_LungMet-1_Rep_3`)) %>%
  dplyr::select(-c(`AKPS_Met-1_Rep_1`, `AKPS_Met-1_Rep_2`, `AKPS_Met-1_Rep_3`)) %>%
  dplyr::select(-c(`AKPS_TdLN-1_Rep_1`, `AKPS_TdLN-1_Rep_2`, `AKPS_TdLN-1_Rep_3`)) %>%
  dplyr::select(-c(`AKP_Tu-2_Rep_1`, `AKP_Tu-2_Rep_2`)) %>% 
  dplyr::select(-c(Normal_002_1, Normal_002_2,Normal_002_3)) 


#####now lets average biological replicates for visualization purposes ONLY which I think is allowed but double check with Justin ###### 
dmmr_pole_tpm_AVERAGED$AK_TU <- (dmmr_pole_tpm_AVERAGED$AK_Tu_1 + dmmr_pole_tpm_AVERAGED$AK_Tu_2 + dmmr_pole_tpm_AVERAGED$AK_Tu_3)/3
dmmr_pole_tpm_AVERAGED$Polyp <- (dmmr_pole_tpm_AVERAGED$Polyp_1 + dmmr_pole_tpm_AVERAGED$Polyps_2 + dmmr_pole_tpm_AVERAGED$Polyps_3)/3
dmmr_pole_tpm_AVERAGED$AKPS_Tu <- (dmmr_pole_tpm_AVERAGED$`AKPS_Tu-1` + dmmr_pole_tpm_AVERAGED$`AKPS_Tu-2` + dmmr_pole_tpm_AVERAGED$`AKPS_Tu-3`)/3
dmmr_pole_tpm_AVERAGED$AKP_Tu <- (dmmr_pole_tpm_AVERAGED$`AKP_Tu-1` + dmmr_pole_tpm_AVERAGED$`AKP_Tu-2`)/2

dmmr_pole_tpm_AVERAGED <- dmmr_pole_tpm_AVERAGED %>% dplyr::select(-c(`AK_Tu_1`,`AK_Tu_2`,`AK_Tu_3`)) %>% 
  dplyr::select(-c(`Polyp_1`,`Polyps_2`,`Polyps_3`)) %>% 
  dplyr::select(-c(`AKPS_Tu-1`,`AKPS_Tu-2`,`AKPS_Tu-3`)) %>%  
  dplyr::select(-c(`AKP_Tu-1`,`AKP_Tu-2`))

####################################stop averaging stuff############################
rownames(dmmr_pole_tpm_AVERAGED) <- NULL
dmmr_pole_tpm_AVERAGED <- dmmr_pole_tpm_AVERAGED %>% column_to_rownames("gene")

res_TPM_normalized <- deconvolute_mouse(dmmr_pole_tpm_AVERAGED, method = 'mmcp_counter')

res_TPM_normalized_long <- res_TPM_normalized %>% pivot_longer(cols=-cell_type, names_to="sample", values_to="count")

res_TPM_normalized_matrix <- res_TPM_normalized_long %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  column_to_rownames(var = "cell_type") %>%
  as.matrix()

colOrder <- c("Normal_002", "Polyp", "AK_TU", "AKP_Tu","AKPS_Tu","AKPS_Met", "AKPS_TdLN", "AKPS_LungMet")
res_TPM_normalized_matrix <- res_TPM_normalized_matrix[,colOrder]

#make heat map 
pdf("heatMap.pdf")
pheatmap(res_TPM_normalized_matrix, 
         cluster_rows = FALSE,  # Cluster rows (cell types)
         cluster_cols = FALSE,  # Cluster columns (samples)
         scale = "row",       # Scale by row (cell type)
         main = "Cell Counts Across Samples",  # Title
         xlab = "Sample",      # X-axis label
         ylab = "Cell Type",    # Y-axis label
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color palette
         fontsize = 10,        # Font size
         fontsize_row = 8,     # Row font size
         fontsize_col = 8, # Column font size
         show_rownames = TRUE
)
dev.off()

#now move on to making line plot
res_TPM_normalized_matrix_long <- as.data.frame(res_TPM_normalized_matrix) %>% rownames_to_column("Cell_Type") %>% 
  pivot_longer(cols=-Cell_Type,names_to = "Stage", values_to = "Fraction")

res_TPM_normalized_matrix_long$Stage <- factor(res_TPM_normalized_matrix_long$Stage, levels = c("Normal_002", "Polyp", "AK_TU", "AKP_Tu", "AKPS_Tu", 
                                                      "AKPS_Met", "AKPS_TdLN", "AKPS_LungMet"))

# Plot
pdf("linePlot.pdf")
ggplot(res_TPM_normalized_matrix_long, aes(x = Stage, y = Fraction, group = Cell_Type, color = Cell_Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Immune Cell Fraction Across Tumor Progression", x = "Stage", y = "Fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

########## create groupings of cells and then move onto making plots ###############


res_TPM_normalized_matrix_long_GROUPED <- res_TPM_normalized_matrix_long %>% mutate(Cell_Type = case_when(
  Cell_Type %in% c("B cell","B cell memory") ~ "B cell",
  Cell_Type %in% c("Monocyte", "Macrophage/Monocyte", "Granulocyte-monocyte progenitor") ~ "Myeloid Cells",
  ! Cell_Type %in% c("B cell","Myeloid Cells") ~ Cell_Type
)) %>% 
  group_by(Stage, Cell_Type) %>% 
  summarise(Fraction = mean(Fraction), .groups="drop")

res_TPM_normalized_matrix_long_GROUPED <- res_TPM_normalized_matrix_long_GROUPED[!res_TPM_normalized_matrix_long_GROUPED$Cell_Type %in% c("Cancer associated fibroblast","Endothelial cell"),]

res_TPM_normalized_matrix_long_GROUPED_summary <- res_TPM_normalized_matrix_long_GROUPED %>%
  group_by(Stage, Cell_Type) %>%
  summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop")

#make lineplot 
pdf("linePlot_grouped.pdf")
ggplot(res_TPM_normalized_matrix_long_GROUPED, aes(x = Stage, y = Fraction, group = Cell_Type, color = Cell_Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Immune Cell Fraction Across Tumor Progression", x = "Stage", y = "Fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#prepare data for heat map 
res_TPM_normalized_matrix_GROUPED <- res_TPM_normalized_matrix_long_GROUPED %>%
  pivot_wider(names_from = Stage, values_from = Fraction, values_fill = 0) %>%
  column_to_rownames("Cell_Type")  # Set Cell_Type as row names

png("heatMap_grouped.png", width=8, height=6, unit="in", res=400)
pheatmap(res_TPM_normalized_matrix_GROUPED, 
         cluster_rows = FALSE,  # Cluster rows (cell types)
         cluster_cols = FALSE,  # Cluster columns (samples)
         scale = "row",       # Scale by row (cell type)
         main = "Immune Cell Counts Across Samples",  # Title
         xlab = "Sample",      # X-axis label
         ylab = "Cell Type",    # Y-axis label
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color palette
         fontsize = 15,        # Font size
         fontsize_row = 14,     # Row font size
         fontsize_col = 14, # Column font size
         show_rownames = TRUE
)
dev.off()

#start making alluvial plot# 
res_TPM_normalized_matrix_long_GROUPED <- res_TPM_normalized_matrix_long_GROUPED %>%
  mutate(Stage = factor(Stage, levels = c(
    "Normal_002", "Polyp", "AK_TU", "AKP_Tu",
    "AKPS_Tu", "AKPS_Met", "AKPS_TdLN", "AKPS_LungMet"
  )))

res_TPM_normalized_matrix_long_GROUPED_summary_subset <- res_TPM_normalized_matrix_long_GROUPED %>%
  group_by(Stage, Cell_Type) %>%
  summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop")

pdf("AlluvialPlot.pdf")

ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_subset,
       aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
           y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cell Type Fractions Across Stages",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

dev.off()

res_TPM_normalized_matrix_long_GROUPED_summary_subset_subset <- res_TPM_normalized_matrix_long_GROUPED_summary_subset[res_TPM_normalized_matrix_long_GROUPED_summary_subset$Cell_Type %in% c("NK cell","Neutrophil","T cell","T cell CD8+"),]

res_TPM_normalized_matrix_long_GROUPED_summary_subset_subset

pdf("AlluvialPlotSubset.pdf")

ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_subset_subset,
       aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
           y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cell Type Fractions Across Stages",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

dev.off()

#make alluvial plots with revisions - start with all immune cells
res_TPM_normalized_matrix_long_GROUPED_summary_1 <- res_TPM_normalized_matrix_long_GROUPED_summary[res_TPM_normalized_matrix_long_GROUPED_summary$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_Met"),]
res_TPM_normalized_matrix_long_GROUPED_summary_2 <- res_TPM_normalized_matrix_long_GROUPED_summary[res_TPM_normalized_matrix_long_GROUPED_summary$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_TdLN"),]
res_TPM_normalized_matrix_long_GROUPED_summary_3 <- res_TPM_normalized_matrix_long_GROUPED_summary[res_TPM_normalized_matrix_long_GROUPED_summary$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_LungMet"),]

p1 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_1,
       aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
           y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

p2 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_2,
       aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
           y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cell Type Fractions Across Stages",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

p3 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_3,
       aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
           y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#make alluvial plots with revisions - move onto subset of immune cells 
res_TPM_normalized_matrix_long_GROUPED_summary_subset_1 <- res_TPM_normalized_matrix_long_GROUPED_summary_subset[res_TPM_normalized_matrix_long_GROUPED_summary_subset$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_Met"),]
res_TPM_normalized_matrix_long_GROUPED_summary_subset_2 <- res_TPM_normalized_matrix_long_GROUPED_summary_subset[res_TPM_normalized_matrix_long_GROUPED_summary_subset$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_TdLN"),]
res_TPM_normalized_matrix_long_GROUPED_summary_subset_3 <- res_TPM_normalized_matrix_long_GROUPED_summary_subset[res_TPM_normalized_matrix_long_GROUPED_summary_subset$Stage %in% c("AK_TU","AKP_Tu","AKPS_Tu","AKPS_LungMet"),]

p1 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_subset_1,
             aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
                 y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

p2 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_subset_2,
             aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
                 y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cell Type Fractions Across Stages",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

p3 <- ggplot(res_TPM_normalized_matrix_long_GROUPED_summary_subset_3,
             aes(x = Stage, stratum = Cell_Type, alluvium = Cell_Type,
                 y = Fraction, fill = Cell_Type)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.7) +
  geom_stratum(alpha = 0.8) +
  theme_minimal() +
  labs(title = "",
       y = "Fraction", x = "Stage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#perform differential gene expression ##### 
library(DESeq2)
dmmr_pole_counts_DDS <- as.data.frame(dmmr_pole_counts)
dmmr_pole_counts_DDS$AKPS_LungMet <- (dmmr_pole_counts_DDS$`AKPS_LungMet-1_Rep_1` +
                                                              dmmr_pole_counts_DDS$`AKPS_LungMet-1_Rep_2` +
                                        dmmr_pole_counts_DDS$`AKPS_LungMet-1_Rep_3`) 

dmmr_pole_counts_DDS$AKPS_Met <- (dmmr_pole_counts_DDS$`AKPS_Met-1_Rep_1` +
                                    dmmr_pole_counts_DDS$`AKPS_Met-1_Rep_2` +
                                    dmmr_pole_counts_DDS$`AKPS_Met-1_Rep_3`) 
dmmr_pole_counts_DDS$AKPS_TdLN <- (dmmr_pole_counts_DDS$`AKPS_TdLN-1_Rep_1` +
                                     dmmr_pole_counts_DDS$`AKPS_TdLN-1_Rep_2` +
                                     dmmr_pole_counts_DDS$`AKPS_TdLN-1_Rep_3`) 

dmmr_pole_counts_DDS$`AKP_Tu-2` <- (dmmr_pole_counts_DDS$`AKP_Tu-2_Rep_1` +
                                      dmmr_pole_counts_DDS$`AKP_Tu-2_Rep_2`)

dmmr_pole_counts_DDS<- dmmr_pole_counts_DDS%>%
  dplyr::select(
    -c(
      `AKPS_Met-1_Rep_1`, `AKPS_Met-1_Rep_2`, `AKPS_Met-1_Rep_3`,
      `AKPS_TdLN-1_Rep_1`, `AKPS_TdLN-1_Rep_2`, `AKPS_TdLN-1_Rep_3`,
      `AKPS_LungMet-1_Rep_1`, `AKPS_LungMet-1_Rep_2`, `AKPS_LungMet-1_Rep_3`,
      `AKP_Tu-2_Rep_1`, `AKP_Tu-2_Rep_2`
    )
  )

#now work on metadata frame 

colonCancerMetaData_GSEA <- data.frame(colnames(dmmr_pole_counts_DDS))

names(colonCancerMetaData_GSEA)[1] <- "sample"
colonCancerMetaData_GSEA$grp <- c(rep("AK_Tu",3), rep("Normal",3), rep("Polyps",3), rep("AKPS_Tu",3), "AKP_Tu", "AKPS_LungMet", "AKPS_Met", "AKPS_TdLN", "AKP_Tu")

dds <- DESeqDataSetFromMatrix(countData = dmmr_pole_counts_DDS,
                              colData = colonCancerMetaData_GSEA,
                              design= ~grp)
dds <- DESeq(dds)

#comparisn of AKvAKP
res_AKvAKP <- lfcShrink(dds,contrast=c("grp","AKP_Tu","AK_Tu"), type='ashr')
res_AKvAKP$ranking <- sign(res_AKvAKP$log2FoldChange) * -log10(res_AKvAKP$pvalue)
res_AKvAKP <- as.data.frame(res_AKvAKP)
res_AKvAKP <- res_AKvAKP %>% rownames_to_column(var = "gene")
res_AKvAKP <- res_AKvAKP[,c("gene","ranking")]
res_AKvAKP <- res_AKvAKP[order(-res_AKvAKP$ranking),]
write.table(res_AKvAKP, file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Preranked Genesets/res_AKvAKP2.rnk", sep = "\t", row.names = FALSE, col.names = c("gene", ""), quote = FALSE)
###comparison of AKPvAKPS 
res_AKPvAKPS <- lfcShrink(dds,contrast=c("grp","AKPS_Tu","AKP_Tu"), type='ashr')
res_AKPvAKPS$ranking <- sign(res_AKPvAKPS$log2FoldChange) * -log10(res_AKPvAKPS$pvalue)
res_AKPvAKPS <- as.data.frame(res_AKPvAKPS)
res_AKPvAKPS <- res_AKPvAKPS %>% rownames_to_column(var = "gene")
res_AKPvAKPS <- res_AKPvAKPS[,c("gene","ranking")]
res_AKPvAKPS <- res_AKPvAKPS[order(-res_AKPvAKPS$ranking),]
write.table(res_AKPvAKPS, file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Preranked Genesets/res_AKPvAKPS.rnk", sep = "\t", row.names = FALSE, col.names = c("gene", ""), quote = FALSE)
###comparison of AKPSvLungMet 
res_AKPSvLungMet <- lfcShrink(dds,contrast=c("grp","AKPS_LungMet","AKPS_Tu"), type='ashr')
res_AKPSvLungMet$ranking <- sign(res_AKPSvLungMet$log2FoldChange) * -log10(res_AKPSvLungMet$pvalue)
res_AKPSvLungMet <- as.data.frame(res_AKPSvLungMet)
res_AKPSvLungMet <- res_AKPSvLungMet %>% rownames_to_column(var = "gene")
res_AKPSvLungMet <- res_AKPSvLungMet[,c("gene","ranking")]
res_AKPSvLungMet <- res_AKPSvLungMet[order(-res_AKPSvLungMet$ranking),]
write.table(res_AKPSvLungMet, file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Preranked Genesets/res_AKPSvLungMet.rnk", sep = "\t", row.names = FALSE, col.names = c("gene", ""), quote = FALSE)
###comparison of AKPSvTDLN 
res_AKPSvTLDN <- lfcShrink(dds, contrast = c("grp", "AKPS_TdLN", "AKPS_Tu"), type = 'ashr')
res_AKPSvTLDN$ranking <- sign(res_AKPSvTLDN$log2FoldChange) * -log10(res_AKPSvTLDN$pvalue)
res_AKPSvTLDN <- as.data.frame(res_AKPSvTLDN)
res_AKPSvTLDN <- res_AKPSvTLDN %>% rownames_to_column(var = "gene")
res_AKPSvTLDN <- res_AKPSvTLDN[, c("gene", "ranking")]
res_AKPSvTLDN <- res_AKPSvTLDN[order(-res_AKPSvTLDN$ranking), ]
write.table(res_AKPSvTLDN, file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Preranked Genesets/res_AKPSvTLDN.rnk", sep = "\t", row.names = FALSE, col.names = c("gene", ""), quote = FALSE)
###comparison of AKPSvMet 
res_AKPSvMet <- lfcShrink(dds, contrast = c("grp", "AKPS_Met", "AKPS_Tu"), type = 'ashr')
res_AKPSvMet$ranking <- sign(res_AKPSvMet$log2FoldChange) * -log10(res_AKPSvMet$pvalue)
res_AKPSvMet <- as.data.frame(res_AKPSvMet)
res_AKPSvMet <- res_AKPSvMet %>% rownames_to_column(var = "gene")
res_AKPSvMet <- res_AKPSvMet[, c("gene", "ranking")]
res_AKPSvMet <- res_AKPSvMet[order(-res_AKPSvMet$ranking), ]
write.table(res_AKPSvMet, file = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Preranked Genesets/res_AKPSvMet.rnk", sep = "\t", row.names = FALSE, col.names = c("gene", ""), quote = FALSE)

######use the dmmr_pole_tpm data frame (TPM normalized counts) to perform GSEA ####### 
dmmr_pole_tpm_GSEA  <- dmmr_pole_tpm %>% rename("gene" = "NAME")
dmmr_pole_tpm_GSEA  <- dmmr_pole_tpm_GSEA %>% add_column(description="na", .after="NAME")
#average technical replicates 
dmmr_pole_tpm_GSEA$AKPS_LungMet <- (dmmr_pole_tpm_GSEA $`AKPS_LungMet-1_Rep_1` +
                                       dmmr_pole_tpm_GSEA $`AKPS_LungMet-1_Rep_2` +
                                          dmmr_pole_tpm_GSEA$`AKPS_LungMet-1_Rep_3`) / 3

dmmr_pole_tpm_GSEA$AKPS_Met <- (dmmr_pole_tpm_GSEA$`AKPS_Met-1_Rep_1` +
                                      dmmr_pole_tpm_GSEA$`AKPS_Met-1_Rep_2` +
                                      dmmr_pole_tpm_GSEA$`AKPS_Met-1_Rep_3`) / 3
dmmr_pole_tpm_GSEA$AKPS_TdLN <- (dmmr_pole_tpm_GSEA$`AKPS_TdLN-1_Rep_1` +
                                       dmmr_pole_tpm_GSEA$`AKPS_TdLN-1_Rep_2` +
                                       dmmr_pole_tpm_GSEA$`AKPS_TdLN-1_Rep_3`) / 3

dmmr_pole_tpm_GSEA$`AKP_Tu-2` <- (dmmr_pole_tpm_GSEA$`AKP_Tu-2_Rep_1` +
                                        dmmr_pole_tpm_GSEA$`AKP_Tu-2_Rep_2`) / 2

dmmr_pole_tpm_GSEA <- dmmr_pole_tpm_GSEA %>%
  dplyr::select(
    -c(
      `AKPS_Met-1_Rep_1`, `AKPS_Met-1_Rep_2`, `AKPS_Met-1_Rep_3`,
      `AKPS_TdLN-1_Rep_1`, `AKPS_TdLN-1_Rep_2`, `AKPS_TdLN-1_Rep_3`,
      `AKPS_LungMet-1_Rep_1`, `AKPS_LungMet-1_Rep_2`, `AKPS_LungMet-1_Rep_3`,
      `AKP_Tu-2_Rep_1`, `AKP_Tu-2_Rep_2`
    )
  )

write.table(dmmr_pole_tpm_GSEA,"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Normalized Counts/normCountsGSEA2.txt", append = TRUE, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

colonCancerMetaData_GSEA <- colonCancerMetaData_GSEA %>% rownames_to_column(var="sample")

##write .cls file 
source("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/CREB5 Breast Cancer/R-Code + Analysis/write_cls_function_for_gsea.R")
write_cls(file_prefix = "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/CLS files",counts = dmmr_pole_tpm_GSEA[,3:19], design=colonCancerMetaData_GSEA, var_to_test = "grp" , libID_col = "sample")


#clean up the results here of AKvAKP
res_AKvAKP_1 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKvAKP.GseaPreranked.1755899800051/gsea_report_for_na_neg_1755899800051.tsv")
res_AKvAKP_2 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKvAKP.GseaPreranked.1755899800051/gsea_report_for_na_pos_1755899800051.tsv") 

res_AKvAKP <- rbind(res_AKvAKP_1, res_AKvAKP_2)
res_AKvAKP <- res_AKvAKP[, c("NAME", "NES", "NOM p-val", "FDR q-val")]
write.xlsx(res_AKvAKP, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Results/res_AKvAKP.xlsx")


#####clean up results here of AKPvAKPS
res_AKPvAKPS_1 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPvAKPS.GseaPreranked.1755899965941/gsea_report_for_na_neg_1755899965941.tsv")
res_AKPvAKPS_2 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPvAKPS.GseaPreranked.1755899965941/gsea_report_for_na_pos_1755899965941.tsv") 

res_AKPvAKPS <- rbind(res_AKPvAKPS_1, res_AKPvAKPS_2)
res_AKPvAKPS <- res_AKPvAKPS[, c("NAME", "NES", "NOM p-val", "FDR q-val")]
write.xlsx(res_AKPvAKPS, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Results/res_AKPvAKPS.xlsx")

#####clean up results here of AKPSvLungMet

res_AKPSvLungMet1 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvLungMet.GseaPreranked.1755900265563/gsea_report_for_na_neg_1755900265563.tsv")
res_AKPSvLungMet2 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvLungMet.GseaPreranked.1755900265563/gsea_report_for_na_pos_1755900265563.tsv") 

res_AKPSvLungMet <- rbind(res_AKPSvLungMet1, res_AKPSvLungMet2)
res_AKPSvLungMet <- res_AKPSvLungMet[, c("NAME", "NES", "NOM p-val", "FDR q-val")]
write.xlsx(res_AKPSvLungMet, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Results/res_AKPSvLungMet.xlsx")


#clean up results of AKPSvTdLN 
res_AKPSvTdLN1 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvTdLN.GseaPreranked.1755900033191/gsea_report_for_na_neg_1755900033191.tsv")
res_AKPSvTdLN2 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvTdLN.GseaPreranked.1755900033191/gsea_report_for_na_pos_1755900033191.tsv") 

res_AKPSvTdLN <- rbind(res_AKPSvTdLN1, res_AKPSvTdLN2)
res_AKPSvTdLN <- res_AKPSvTdLN[, c("NAME", "NES", "NOM p-val", "FDR q-val")]
write.xlsx(res_AKPSvTdLN, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Results/res_AKPSvTdLN.xlsx")

###clean up results of AKPSvMet 
res_AKPSvMet1 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvMet.GseaPreranked.1755900147571/gsea_report_for_na_neg_1755900147571.tsv")
res_AKPSvMet2 <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Output/AKPSvMet.GseaPreranked.1755900147571/gsea_report_for_na_pos_1755900147571.tsv") 

res_AKPSvMet <- rbind(res_AKPSvMet1, res_AKPSvMet2)
res_AKPSvMet <- res_AKPSvMet[, c("NAME", "NES", "NOM p-val", "FDR q-val")]
write.xlsx(res_AKPSvMet, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/GSEA/Results/res_AKPSvMet.xlsx")


###make heatmap showing the different pathways ###### 
setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output")
prepare_heatmap_data <- function(df, comparison_name) {
  df %>%
    dplyr::select(NAME, NES, `FDR q-val`) %>%
    mutate(Comparison = comparison_name,
           Significance = case_when(
             `FDR q-val` <= 0.05 ~ "*",
             `FDR q-val` <= 0.1 ~ "+",
             TRUE ~ ""
           ))
}

ak_akp_data <- prepare_heatmap_data(res_AKvAKP, "AK vs. AKP")
akp_akps_data <- prepare_heatmap_data(res_AKPvAKPS, "AKP vs. AKPS")
akps_TdLN_data <- prepare_heatmap_data(res_AKPSvTdLN, "AKPS vs. AKPS_TdLN")
akps_LungMet_data <- prepare_heatmap_data(res_AKPSvLungMet, "AKPS vs. AKPS_LungMet")
akps_DistMet_data <- prepare_heatmap_data(res_AKPSvMet, "AKPS vs. AKPS_DistMet")

combined_data <- bind_rows(ak_akp_data, akp_akps_data, akps_TdLN_data, akps_LungMet_data,akps_DistMet_data)

# 3. Categorize Pathways and Order

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

combined_data <- combined_data %>%
  mutate(
    Pathway_Type = ifelse(NAME %in% immune_pathways, "Immune", "Non-Immune"),
    # Create a numeric index for ordering within each category
    Pathway_Order = ifelse(Pathway_Type == "Immune",
                           match(NAME, immune_pathways),
                           match(NAME, setdiff(unique(NAME), immune_pathways)))
  ) %>%
  arrange(desc(Pathway_Type), Pathway_Order) %>%  # Changed this line
  mutate(NAME = factor(NAME, levels = unique(NAME))) # Crucial for ordering

# 4. Create the Heatmaps with ggplot2 and Patchwork
# Common heatmap aesthetics

combined_data_immune <- combined_data[combined_data$NAME %in% immune_pathways,]

base_heatmap_immune <- ggplot(combined_data_immune, aes(Comparison, NAME, fill = NES)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 8),
    strip.text = element_text(size = 11, face = "bold")
  ) +
  facet_wrap(~Comparison, nrow = 1, scales = "free_x")

# Add significance stars
final_heatmap_immune <- base_heatmap_immune +
  geom_text(aes(label = Significance), color = "black", size = 9) +
  theme(
    strip.text.x = element_text(size = 11, margin = ggplot2::margin(t = 3, b = 3)),  # Bigger text, more margin
  )

# Save the combined plot
ggsave("combined_pathway_heatmap_FINAL_immune.pdf", final_heatmap_immune, width = 14, height = 8)


### repeat same code as immune pathways but for oncogenic pathways #### 
combined_data_oncogenic<- combined_data[combined_data$NAME %in% oncogenic_pathways,]

base_heatmap_oncogenic <- ggplot(combined_data_oncogenic, aes(Comparison, NAME, fill = NES)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 8), 
    strip.text = element_text(size = 9, face = "bold")
  ) 

# Add significance stars
final_heatmap_oncogenic <- base_heatmap_oncogenic +
  geom_text(aes(label = Significance), color = "black", size = 9) +
  theme(
    strip.text.x = element_text(size = 9, margin = ggplot2::margin(t = 3, b = 3)),  # Bigger text, more margin
  )

# Save the combined plot
ggsave("combined_pathway_heatmap_FINAL_oncogenic.pdf", final_heatmap_oncogenic, width = 14, height = 8)



######now use dmmr_pole_tpm to utilize the MmCMS package to assign CMS subtypes to mouse colorectal cancer samples
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install MmCMS package
devtools::install_github("MolecularPathologyLab/MmCMS")

# Load the package
library(MmCMS)

dmmr_pole_tpm_CMS <- dmmr_pole_tpm %>% column_to_rownames("gene")
dmmr_pole_tpm_CMS <-log2(dmmr_pole_tpm_CMS + 1) 
#check with Ella to see if taking log of TPM before running MmCMS is appropriate 

#now some shenangians to limit us to MGI gene symbols to protect the p values of MmCMS classifier 

mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror="useast")
valid_mgi_symbols <- getBM(
  attributes = c("mgi_symbol"),
  mart = mouse
) %>%
  filter(mgi_symbol != "") %>%
  pull(mgi_symbol) %>%
  unique()

### stop those shenanignas and start to filter 
dmmr_pole_tpm_CMS_clean <- dmmr_pole_tpm_CMS[rownames(dmmr_pole_tpm_CMS) %in% valid_mgi_symbols, ]

results_CMS_types <- MmCMS(as.matrix(dmmr_pole_tpm_CMS_clean), templates = MmCMS::template.CMS.C, Genesets=c("template.CMS.C"), seed=1, FDR=1)


library(dplyr)
library(tidyr)
library(ggplot2)

# Assuming your data is in a data frame called `cms_data`
# And the sample names are in rownames
cms_data <- results_CMS_types # Replace with your actual data object

# Add sample names as a column
cms_data <- cms_data %>%
  rownames_to_column("sample")

# Extract the group from the sample name (e.g., AK_Tu_1 → AK_Tu)
cms_data <- cms_data %>%
  mutate(group = sub("(_\\d+.*|-[^_]+)$", "", sample))  # Strips replicate or sample suffixes

cms_data$group <- ifelse(cms_data$group == "AKP_Tu-2_Rep", "AKP_Tu",cms_data$group)
cms_data$group <- ifelse(cms_data$group == "AKPS_LungMet-1_Rep", "AKPS LungMet",cms_data$group)
cms_data$group <- ifelse(cms_data$group == "AKPS_Met-1_Rep", "AKPS Met",cms_data$group)
cms_data$group <- ifelse(cms_data$group == "AKPS_TdLN-1_Rep", "AKPS TdLN",cms_data$group)
cms_data$group <- ifelse(cms_data$group == "Polyps", "Polyp",cms_data$group)
#cms_data$prediction <- ifelse(is.na(cms_data$prediction), "Prediction Not Stat. Sig.",cms_data$prediction)
cms_data$prediction <- ifelse(is.na(cms_data$prediction), "Unclassified", cms_data$prediction)
cms_counts <- cms_data %>%
  dplyr::count(group, prediction)

group_order <- c("Normal", "Polyp", "AK_Tu", "AKP_Tu", "AKPS_Tu", "AKPS Met", "AKPS LungMet", "AKPS TdLN")

# Set group as a factor with specified order
cms_counts$group <- factor(cms_counts$group, levels = group_order)
cms_counts$prediction <-factor(cms_counts$prediction) 
stackedBar_Mouse_CMS <- ggplot(cms_counts, aes(x = group, y = n, fill = prediction)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1", na.value = "black") +
  labs(
    title = "CMS Subtype Distribution by Group (Counts)",
    x = "Group",
    y = "Number of Samples",
    fill = "CMS Subtype"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/stackedBar_Mouse_CMS.pdf",stackedBar_Mouse_CMS, width = 10, height = 8)


###make a volcano plot here for figure 2C #### 

get_top_genes <- function(df, fc_filter, n = 10) {
  subset_df <- df[fc_filter(df$log2FoldChange), ]
  subset_df <- subset_df[order(subset_df$fdr_adj_pvalue), ]
  rownames(subset_df)[1:min(n, nrow(subset_df))]
}

# Function: shrinkage, adj p-values, top genes, and return ggplot object
make_volcano <- function(dds, contrast, plot_title, extra_labels = NULL,
                         fc_cut = 2, p_cut = 0.05, n_label = 10,
                         show_legend = TRUE) {
  
  res <- lfcShrink(dds, contrast = contrast, type = 'ashr')
  res$fdr_adj_pvalue <- p.adjust(res$pvalue, method = "fdr")
  
  up_genes <- get_top_genes(res, function(fc) fc >= fc_cut, n = n_label)
  down_genes <- get_top_genes(res, function(fc) fc <= -fc_cut, n = n_label)
  
  select_labels <- unique(c(up_genes, down_genes, extra_labels, "Reg4"))
  
  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    selectLab = select_labels,
    x = 'log2FoldChange',
    y = 'fdr_adj_pvalue',
    title = plot_title,
    pCutoff = p_cut,
    FCcutoff = fc_cut,
    drawConnectors = TRUE,
    maxoverlapsConnectors = Inf,
    lengthConnectors = unit(0.0075, "npc"),
    boxedLabels = FALSE
  )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# Create plots
p1 <- make_volcano(dds, c("grp", "AKP_Tu", "AK_Tu"), "AK vs AKP", show_legend = TRUE)
p2 <- make_volcano(dds, c("grp", "AKPS_Tu", "AKP_Tu"), "AKP vs AKPS", show_legend = FALSE)
p3 <- make_volcano(dds, c("grp", "AKPS_LungMet", "AKPS_Tu"), "AKPS vs LungMet", show_legend = FALSE)
p4 <- make_volcano(dds, c("grp", "AKPS_TdLN", "AKPS_Tu"), "AKPS vs TdLN", show_legend = FALSE)
p5 <- make_volcano(dds, c("grp", "AKPS_Met", "AKPS_Tu"), "AKPS vs DistMet", show_legend = FALSE)

# Extract legend from first plot
legend <- get_legend(p1 + theme(legend.position = "bottom"))

# Remove legend from p1 for plotting
p1_clean <- p1 + theme(legend.position = "none")

# Arrange plots (3x2 layout) with single legend below
pdf("Volcano_Plots_AllComparisons_Dec12.pdf", width = 16, height = 12)
combined_plots <- plot_grid(p1_clean, p2, p3, p4, p5, ncol = 3, align = "hv")
final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.08))
print(final_plot)
dev.off()

################## now do some work to prepare the various significant genes as a starting point for cluster analysis 
View(assay(dds))
View(colonCancerMetaData_GSEA)

##########get list of genes unique to AK ############# 

rownames(colonCancerMetaData_GSEA) <- colonCancerMetaData_GSEA$sample

dmmr_pole_counts_DDS <- dmmr_pole_counts_DDS[, 
                                             (colnames(dmmr_pole_counts_DDS) %in% c("AK_Tu_1", "AK_Tu_2","AK_Tu_3", "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3", "AKP_Tu-1", "AKP_Tu-2")),]

colonCancerMetaData_GSEA <- colonCancerMetaData_GSEA[
  (colonCancerMetaData_GSEA$sample %in% c("AK_Tu_1", "AK_Tu_2","AK_Tu_3", "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3", "AKP_Tu-1", "AKP_Tu-2") ),
]

colonCancerMetaData_GSEA$grp_relevel_AK <- ifelse(colonCancerMetaData_GSEA$grp == "AK_Tu", "AK_Tu", "Other")
colonCancerMetaData_GSEA$grp_relevel_AK <- factor(colonCancerMetaData_GSEA$grp_relevel_AK, levels = c("Other","AK_Tu"))

dds_AK <- DESeqDataSetFromMatrix(countData = dmmr_pole_counts_DDS,
                              colData = colonCancerMetaData_GSEA,
                              design= ~grp_relevel_AK)
dds_AK <- DESeq(dds_AK)


res_AK <- lfcShrink(
  dds_AK,
  contrast = c("grp_relevel_AK", "AK_Tu", "Other"),
  type = "ashr"
)


############ get list of genes unique to AKP ############## 

colonCancerMetaData_GSEA$grp_relevel_AKP <- ifelse(colonCancerMetaData_GSEA$grp == "AKP_Tu", "AKP_Tu", "Other")
colonCancerMetaData_GSEA$grp_relevel_AKP <- factor(colonCancerMetaData_GSEA$grp_relevel_AKP, levels = c("Other","AKP_Tu"))

dds_AKP <- DESeqDataSetFromMatrix(countData = dmmr_pole_counts_DDS,
                                 colData = colonCancerMetaData_GSEA,
                                 design= ~grp_relevel_AKP)
dds_AKP <- DESeq(dds_AKP)


res_AKP <- lfcShrink(
  dds_AKP,
  contrast = c("grp_relevel_AKP", "AKP_Tu", "Other"),
  type = "ashr"
)


############ get list of genes unique to AKPS ############## 
colonCancerMetaData_GSEA$grp_relevel_AKPS <- ifelse(colonCancerMetaData_GSEA$grp == "AKPS_Tu", "AKPS_Tu", "Other")
colonCancerMetaData_GSEA$grp_relevel_AKPS <- factor(colonCancerMetaData_GSEA$grp_relevel_AKPS, levels = c("Other","AKPS_Tu"))

dds_AKPS <- DESeqDataSetFromMatrix(countData = dmmr_pole_counts_DDS,
                                  colData = colonCancerMetaData_GSEA,
                                  design= ~grp_relevel_AKPS)
dds_AKPS <- DESeq(dds_AKPS)


res_AKPS <- lfcShrink(
  dds_AKPS,
  contrast = c("grp_relevel_AKPS", "AKPS_Tu", "Other"),
  type = "ashr"
)

#########
res_AK <- as.data.frame(res_AK)
res_AK <- res_AK[res_AK$padj <=0.05, ]
res_AK <- res_AK[complete.cases(res_AK), ]

res_AK <- res_AK[order(abs(res_AK$log2FoldChange), decreasing = TRUE),] 

AK_genes <- rownames(res_AK)

###########

res_AKP <- as.data.frame(res_AKP)
res_AKP <- res_AKP[res_AKP$padj <=0.05, ]
res_AKP <- res_AKP[complete.cases(res_AKP), ]

res_AKP <- res_AKP[order(abs(res_AKP$log2FoldChange), decreasing = TRUE),] 

AKP_genes <- rownames(res_AKP)

###########

res_AKPS <- as.data.frame(res_AKPS)
res_AKPS <- res_AKPS[res_AKPS$padj <=0.05, ]
res_AKPS <- res_AKPS[complete.cases(res_AKPS), ]

res_AKPS <- res_AKPS[order(abs(res_AKPS$log2FoldChange), decreasing = TRUE),] 

AKPS_genes <- rownames(res_AKPS)


top_AKPS <- AKPS_genes

# Remove AK genes from AKP list
top_AKP <- setdiff(AKP_genes, top_AKPS)

# Remove AK and AKP genes from AKPS list
top_AK <- setdiff(AK_genes, c(top_AKP, top_AKPS))

# Now take top 100 from each (if enough remain)
top_AK <- top_AK[1:min(500, length(top_AK))]
top_AKP <- top_AKP[1:min(500, length(top_AKP))]
top_AKPS <- top_AKPS[1:min(500, length(top_AKPS))]

top_genes <- c(top_AK, top_AKP, top_AKPS)

write.xlsx(as.data.frame(top_genes), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/top_genes.xlsx")


#### now prepare to get data frames of 500 x 2 for each group so we can make our correlation matrix ####

res_AK_small <- res_AK %>% head(1000) %>% dplyr::select(log2FoldChange)
res_AKP_small <- res_AKP %>% head(1000) %>% dplyr::select(log2FoldChange)
res_AKPS_small <- res_AKPS %>% head(1000) %>% dplyr::select(log2FoldChange)


write.xlsx(as.data.frame(res_AK_small), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AK_small.xlsx", rowNames=TRUE)
write.xlsx(as.data.frame(res_AKP_small), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AKP_small.xlsx",rowNames=TRUE)
write.xlsx(as.data.frame(res_AKPS_small), "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/res_AKPS_small.xlsx",rowNames=TRUE)




##### prepare data for cross-species classifier (on human file)
# FULL PIPELINE: Mouse → Human Cross-species Classification
########################################################

########################################################
# 1. Define samples and metadata
########################################################
comparingSamples <- c(
  "AK_Tu_1", "AK_Tu_2", "AK_Tu_3",       # AK
  "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3", # AKPS
  "AKP_Tu-1", "AKP_Tu-2_Rep_1", "AKP_Tu-2_Rep_2" # AKP
)

CrossSpecComp <- dmmr_pole_counts[, comparingSamples]

bio_ids <- c(
  "AK_Tu-1", "AK_Tu-2", "AK_Tu-3",       # AK
  "AKPS_Tu-1", "AKPS_Tu-2", "AKPS_Tu-3", # AKPS
  "AKP_Tu-1", "AKP_Tu-2", "AKP_Tu-2"     # AKP (rep)
)

groups <- c(rep("AK", 3), rep("AKPS", 3), rep("AKP", 3))

sample_metadata <- data.frame(
  sample_id = colnames(CrossSpecComp),
  bio_id = bio_ids,
  group = groups
)
rownames(sample_metadata) <- colnames(CrossSpecComp)


########################################################
# 2. DESeq2 normalization
########################################################
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = CrossSpecComp,
  colData = sample_metadata,
  design = ~ group
)

dds_collapsed <- collapseReplicates(
  dds,
  groupby = dds$bio_id,
  run = dds$sample_id
)

dds_collapsed <- estimateSizeFactors(dds_collapsed)
vst_CrossSpecComp <- vst(dds_collapsed, blind = TRUE)
CrossSpecComp_norm <- assay(vst_CrossSpecComp)  # genes x samples


########################################################
# 3. Select 8,000 most variable genes
########################################################
gene_vars <- apply(CrossSpecComp_norm, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:8000]
expr_top <- CrossSpecComp_norm[top_genes, ]


########################################################
# 4. Hierarchical clustering (Spearman + Ward.D2)
########################################################
library(dynamicTreeCut)
library(stats)

# Compute correlation distance
cor_mat <- cor(t(expr_top), method="spearman")
dist_mat <- as.dist(1 - cor_mat)

# Hierarchical clustering
hc <- hclust(dist_mat, method = "ward.D2")

# Cut tree dynamically
clusters <- cutreeDynamic(
  dendro = hc,
  distM = as.matrix(dist_mat),
  minClusterSize = 60, 
  method="hybrid"
)

names(clusters) <- hc$labels

# Organize clusters
gene_clusters <- split(names(clusters), clusters)
gene_clusters <- gene_clusters[names(gene_clusters) != "0"]


########################################################
# 5. Filter modules by within-cluster correlation > 0.6
########################################################
keep_modules <- list()

for (m in names(gene_clusters)) {
  gset <- gene_clusters[[m]]
  if (length(gset) < 60) next
  
  sub_expr <- expr_top[gset, ]
  cor_sub <- cor(t(sub_expr), method = "pearson")
  avg_cor <- mean(cor_sub[upper.tri(cor_sub)])
  
  if (avg_cor > 0.6) {
    keep_modules[[paste0("Module_", m)]] <- gset
  }
}

cat(length(keep_modules), "modules retained with r > 0.6 and ≥ 60 genes\n")




#################################################
# 5B. Calculate human versions of modules

human_gene_modules <- list()

for (m in names(gene_clusters)) {
  mouse_genes <- gene_clusters[[m]]
  
  # Map mouse genes to human
  mapped <- map_mouse_to_human_expr(mouse_genes)
  human_genes <- mapped$human_entrez  # extract human gene IDs
  
  # Skip modules with <50% successful mapping
  if (length(human_genes) / length(mouse_genes) < 0.5) next
  
  # Add to human modules list
  human_gene_modules[[paste0("Module_", m)]] <- human_genes
}

############################################################


########################################################
# 6. Compute eigengenes for retained modules
########################################################
library(WGCNA)

# Assuming expr_top is mouse expression (rows = mouse genes)
mapped_expr <- map_mouse_to_human_expr(rownames(expr_top))

# Keep only successfully mapped genes
expr_top_human <- expr_top[rownames(expr_top) %in% mapped_expr$mouse_symbol, ]

# Rename rows to human Entrez IDs
rownames(expr_top_human) <- mapped_expr$human_entrez[match(rownames(expr_top_human), mapped_expr$mouse_symbol)]

module_eigengenes <- list()
for (m in names(human_gene_modules)) {
  gset <- intersect(human_gene_modules[[m]], rownames(expr_top_human))
  if (length(gset) >= 10) {
    module_data <- expr_top_human[gset, , drop = FALSE]
    eig <- prcomp(t(module_data), scale. = TRUE)$x[, 1]
    module_eigengenes[[m]] <- eig
  }
}

# Convert to matrix for downstream use
module_eigengenes <- as.data.frame(module_eigengenes)


########################################################
# 7. Map mouse genes to human orthologs
########################################################
# Assume you already have:
# map_mouse_to_human_expr() → returns data.frame(mouse_gene, human_entrez)

human_gene_modules <- list()

for (m in names(keep_modules)) {
  mouse_genes <- keep_modules[[m]]
  mapped <- map_mouse_to_human_expr(mouse_genes)
  mapped <- mapped[!is.na(mapped$human_entrez), ]
  
  mapping_rate <- nrow(mapped) / length(mouse_genes)
  if (mapping_rate >= 0.5) {
    human_gene_modules[[m]] <- mapped$human_entrez
  }
}

cat(length(human_gene_modules), "modules retained after 50% mapping filter\n")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("GO.db")
#BiocManager::install("impute")

library(WGCNA)

########################################################
# 8. Compute eigengenes for human data
########################################################
# Assuming human_expr (genes x samples) with gene IDs matching human_entrez
TCGA_zscores <- readr::read_delim("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/coad_tcga_gdc/data_mrna_seq_fpkm_zscores_ref_all_samples.txt", delim = "\t") 
TCGA_zscores <- TCGA_zscores %>%
  filter(!is.na(Entrez_Gene_Id)) %>%
  group_by(Entrez_Gene_Id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Entrez_Gene_Id")


TCGA_zscores <- as.matrix(TCGA_zscores)
TCGA_zscores <- TCGA_zscores[complete.cases(TCGA_zscores), ]

eigengenes_human_mat <- sapply(human_gene_modules, function(gset) {
  # keep only genes present in TCGA data
  gset <- intersect(gset, rownames(TCGA_zscores))
  if (length(gset) < 10) return(rep(NA, ncol(TCGA_zscores)))  # skip small modules
  
  sub_expr <- TCGA_zscores[gset, , drop = FALSE]
  # one color per gene so we can compute a single module eigengene
  eig <- WGCNA::moduleEigengenes(t(sub_expr), colors = rep(1, nrow(sub_expr)))$eigengenes[, 1]
  return(eig)
})

eigengenes_human_mat <- as.data.frame(eigengenes_human_mat)
rownames(eigengenes_human_mat) <- colnames(TCGA_zscores)  # samples

######################
#9. Train model on mouse genes 
########################################
mouse_train <- as.data.frame(module_eigengenes)
mouse_train$group <- c("AK", "AK", "AK", "AKP", "AKP", "AKPS", "AKPS", "AKPS")

# Scale predictors
mouse_train_scaled <- mouse_train
mouse_train_scaled[, !names(mouse_train_scaled) %in% "group"] <-
  scale(mouse_train_scaled[, !names(mouse_train_scaled) %in% "group"])

library(nnet)
multi_mod <- multinom(group ~ ., data = mouse_train_scaled)
summary(multi_mod)

######################
#10. Predict human sample labels from mouse model 
########################################

common_mods <- intersect(
  colnames(mouse_train)[-ncol(mouse_train)],
  colnames(eigengenes_human_mat)
)

human_df_scaled <- eigengenes_human_mat[, common_mods, drop = FALSE]
human_df_scaled <- as.data.frame(scale(human_df_scaled))

# 1. Identify all modules the model expects
model_mods <- colnames(mouse_train_scaled)[!colnames(mouse_train_scaled) %in% "group"]

# 2. Subset and reorder human matrix to match model features
missing_mods <- setdiff(model_mods, colnames(eigengenes_human_mat))
human_df_scaled <- eigengenes_human_mat[, intersect(model_mods, colnames(eigengenes_human_mat)), drop = FALSE]

# 3. Add back missing columns as zeros
for (m in missing_mods) {
  human_df_scaled[[m]] <- 0
}

# 4. Reorder columns to match model exactly
human_df_scaled <- human_df_scaled[, model_mods, drop = FALSE]

# 5. scale columns safely -- that is to say never /0

scale_safe <- function(x) {
  if (sd(x, na.rm = TRUE) == 0 || all(is.na(x))) return(rep(0, length(x)))
  as.numeric(scale(x))
}

human_df_scaled <- as.data.frame(apply(human_df_scaled, 2, scale_safe))
rownames(human_df_scaled) <- rownames(eigengenes_human_mat)

# 7. Predict
human_pred_probs <- predict(multi_mod, newdata = human_df_scaled, type = "probs")
human_pred_labels <- predict(multi_mod, newdata = human_df_scaled, type = "class")

# 8. Combine results
pred_results <- data.frame(
  Sample = rownames(human_df_scaled),
  Predicted_Class = human_pred_labels,
  human_pred_probs
)

print(pred_results)


pred_results$highest_pred <- apply(pred_results[,3:5], MARGIN=1, FUN = max)

pred_results$Predicted_Class <- as.factor(pred_results$Predicted_Class)


pred_results$Predicted_Class <- ifelse(pred_results$highest_pred < 0.7, "None", pred_results$Predicted_Class)



pred_results$Predicted_Class <- dplyr::recode(
  pred_results$Predicted_Class,
  "1" = "AK",
  "2" = "AKP",
  "3" = "AKPS"
)

# Keep 'None' as-is
pred_results$Predicted_Class[is.na(pred_results$Predicted_Class)] <- "None"

write.csv(pred_results,"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/all_pred_results.csv")


pred_results <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/all_pred_results.csv")

###write csv files ###
write.csv(pred_results[pred_results$Predicted_Class == "AK",],"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AK_transcriptomic.csv")
write.csv(pred_results[pred_results$Predicted_Class == "AKP",],"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKP_transcriptomic.csv")
write.csv(pred_results[pred_results$Predicted_Class == "AKPS",],"/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKPS_transcriptomic.csv")



##### now make PCA plot to see if these new labels result in better seperation of samples#### 

library(ggplot2)

# ------------------------------------------------------------
# 1. Join predicted class to TCGA expression data
# ------------------------------------------------------------

# Ensure the rownames of TCGA_zscores correspond to genes
# and colnames correspond to samples (as you described)
# Transpose so that samples are rows (for PCA)
TCGA_zscores_t <- t(TCGA_zscores)

#make Sample column 

TCGA_zscores_t <- as.data.frame(TCGA_zscores_t) %>% rownames_to_column("Sample")

pred_results$Sample <- substr(pred_results$Sample,1,15)

# Check matching sample names
common_samples <- intersect(TCGA_zscores_t$Sample, pred_results$Sample)

TCGA_zscores_t <- TCGA_zscores_t[TCGA_zscores_t$Sample %in% common_samples,]

pca_input <- left_join(TCGA_zscores_t, pred_results[, c("Sample", "Predicted_Class")], by = "Sample")

# ------------------------------------------------------------
# 2. PCA on z-score data
# ------------------------------------------------------------

numeric_data <- pca_input[, !names(pca_input) %in% c("Predicted_Class", "Sample")]

keep_cols <- sapply(numeric_data, function(x) sd(x, na.rm = TRUE) > 0)

numeric_data <- numeric_data[, which(keep_cols), drop = FALSE]

rownames(numeric_data) <- TCGA_zscores_t$Sample

pca_res <- prcomp(numeric_data, center = TRUE, scale. = TRUE)

library(ggplot2)

# --- 1. Prepare PCA data frame ---
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

# --- 2. Merge with predicted classes ---
# This ensures samples and predicted classes are aligned correctly
pca_df <- merge(
  pca_df,
  pred_results[, c("Sample", "Predicted_Class")],
  by = "Sample",
  all.x = TRUE
)

# --- 3. Define color mapping ---
color_map <- c(
  "AK" = "red",
  "AKP" = "green",
  "AKPS" = "blue",
  "None" = "grey"
)

# --- 4. Plot PCA with ggplot2 ---
ggplot(pca_df, aes(x = PC1, y = PC2, color = Predicted_Class)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = color_map, na.value = "black") +
  theme_minimal(base_size = 14) +
  xlim(0,50) + 
  ylim(-50,50) +
  labs(
    title = "PCA of TCGA Samples Colored by Predicted Class",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "% variance)")
  ) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

############## create UMAP ######################## 

if (!require(umap)) install.packages("umap")
library(umap)
library(ggplot2)

# numeric_data: rows = samples, columns = genes
# Make sure rows are samples
numeric_data <- t(numeric_data)  # now rows = samples, columns = features

# Run UMAP
set.seed(123)  # for reproducibility
umap_res <- umap(numeric_data)

# Prepare data frame for plotting
umap_df <- data.frame(
  Sample = rownames(numeric_data),
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2]
)


umap_df <- merge(
  umap_df,
  pred_results[, c("Sample", "Predicted_Class")],
  by = "Sample",
  all.x = TRUE
)

# Join predicted classes
umap_df$Predicted_Class <- factor(umap_df$Predicted_Class,
                                  levels = c("AK", "AKP", "AKPS", "None"))

# Color palette
class_colors <- c(
  "AK"   = "#E41A1C",
  "AKP"  = "#377EB8",
  "AKPS" = "#4DAF4A",
  "None" = "grey70"
)

UMAP_Transcriptomic <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Predicted_Class)) +
  geom_point(size = 3, alpha = 0.8) +
  
  # 95% CI ellipse, computed per class
  stat_ellipse(aes(fill = Predicted_Class),
               type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.15, color = NA) +
  
  scale_color_manual(values = class_colors) +
  scale_fill_manual(values = class_colors) +
  
  theme_minimal() +
  labs(
    title = "",
    x = "UMAP1",
    y = "UMAP2",
    color = "Mutation Group",
    fill = "Mutation Group"
  ) +
  # Larger legend
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.5, "lines")
  )

UMAP_Transcriptomic

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/")
ggsave("umap_transcriptomic_Nov9.pdf", UMAP_Transcriptomic)


######### now make a heatmap for quality control ############## 
library(dplyr)
library(pheatmap)

# Merge predicted classes with eigengenes


eigengenes_annot <- eigengenes_human_mat %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(pred_results, by = c("Sample" = "Sample")) %>%
  mutate(Predicted_Class = factor(Predicted_Class,
                                  levels = c("AK", "AKP", "AKPS", "None"))) %>%
  arrange(Predicted_Class)

#remove the ones that don't have a classification ####### 

eigengenes_annot <- eigengenes_annot[eigengenes_annot$Predicted_Class != "None",]


# Extract only module columns (Module1, Module2, etc.)
mat_ordered <- eigengenes_annot %>%
  dplyr::select(starts_with("Module")) %>%
  as.matrix()

rownames(mat_ordered) <- eigengenes_annot$Sample

# Create annotation for row colors
row_anno <- data.frame(Predicted_Class = eigengenes_annot$Predicted_Class)
rownames(row_anno) <- eigengenes_annot$Sample

# Define color palette
ann_colors <- list(
  Predicted_Class = c(AK = "#1b9e77", AKP = "#d95f02",
                      AKPS = "#7570b3")
)

# Make the heatmap
pheatmap(
  mat_ordered,
  cluster_rows = FALSE,       # keep your order (AK → AKP → AKPS → none)
  cluster_cols = TRUE,        # cluster modules if you want
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  scale = "row",              # makes expression patterns easier to compare
  show_rownames = FALSE,
  border_color = NA,
  fontsize_col = 11,
  main = "Module Expression by Predicted Class",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)


######### now move onto making the confusion matrix ############## 
AK_TCGA_genomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AK.csv")
AKP_TCGA_genomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKP.csv")
AKPS_TCGA_genomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/AKPS.csv")

AK_samples_genomic <- AK_TCGA_genomic$Sample.ID
AKP_samples_genomic <- AKP_TCGA_genomic$Sample.ID
AKPS_samples_genomic <- AKPS_TCGA_genomic$Sample.ID

AK_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AK_transcriptomic.csv")
AK_samples_transcriptomic <- str_sub(AK_samples_transcriptomic$Sample, 1, -2)

AKP_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKP_transcriptomic.csv")
AKP_samples_transcriptomic <- str_sub(AKP_samples_transcriptomic$Sample, 1, -2)

AKPS_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKPS_transcriptomic.csv")
AKPS_samples_transcriptomic <- str_sub(AKPS_samples_transcriptomic$Sample, 1, -2)

None_samples_transcriptomic <- pred_results[pred_results$Predicted_Class == "None", "Sample"]

None_samples_transcriptomic <- str_sub(None_samples_transcriptomic, 1, -2)

genomic_df <- data.frame(
  Sample = c(AK_samples_genomic, AKP_samples_genomic, AKPS_samples_genomic),
  Genomic_Class = c(
    rep("AK", length(AK_samples_genomic)),
    rep("AKP", length(AKP_samples_genomic)),
    rep("AKPS", length(AKPS_samples_genomic))
  )
)

# Combine transcriptomic classifications
transcriptomic_df <- data.frame(
  Sample = c(AK_samples_transcriptomic, AKP_samples_transcriptomic, AKPS_samples_transcriptomic, None_samples_transcriptomic),
  Transcriptomic_Class = c(
    rep("AK", length(AK_samples_transcriptomic)),
    rep("AKP", length(AKP_samples_transcriptomic)),
    rep("AKPS", length(AKPS_samples_transcriptomic)), 
    rep("None", length(None_samples_transcriptomic))
    
  )
)

# Merge into one table
class_compare <- merge(genomic_df, transcriptomic_df, by = "Sample", all = TRUE)

class_compare <- class_compare %>%
  filter(Sample %in% genomic_df$Sample)

table_compare <- table(class_compare$Genomic_Class, class_compare$Transcriptomic_Class)

table_compare <- table_compare[, c("None", "AK", "AKP", "AKPS")]

png("ConfusionMatrix_Nov6.png", width = 5, height = 5, units="in", res=300)

pheatmap(
  table_compare,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("white", "steelblue"))(50),
  fontsize = 17,
  fontsize_number = 15,
  angle_col = 0, 
  number_format = "%.0f",
  cellwidth = 60,   # smaller width
  cellheight = 60 
)

dev.off()

# Activate the heatmap viewport
grid.newpage() # start a fresh page
pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.draw(ph$gtable)  # draw the pheatmap object

# Add x-axis label (below the plot)
grid.text("Transcriptomic Classification",
          x = 0.5, y = unit(0.02, "npc"), gp = gpar(fontsize = 13))

# Add y-axis label (rotated, left of the plot)
grid.text("Genomic Classification",
          x = unit(0.02, "npc"), y = 0.5, rot = 90, gp = gpar(fontsize = 13))


############### perform GSEA analysis with the transcriptomic classifications ############# 
AK_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AK_transcriptomic.csv")
AK_samples_transcriptomic <- str_sub(AK_samples_transcriptomic$Sample, 1, -2)

AKP_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKP_transcriptomic.csv")
AKP_samples_transcriptomic <- str_sub(AKP_samples_transcriptomic$Sample, 1, -2)

AKPS_samples_transcriptomic <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/AKPS_transcriptomic.csv")
AKPS_samples_transcriptomic <- str_sub(AKPS_samples_transcriptomic$Sample, 1, -2)


########## now prepare COAD data ##############

message("1. Downloading TCGA-COAD data...")
COAD_query <- GDCquery(project = 'TCGA-COAD',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       workflow.type = 'STAR - Counts')

GDCdownload(query = COAD_query)

COAD_data_obj <- GDCprepare(query = COAD_query, summarizedExperiment = TRUE)

message("2. Extracting gene annotations...")
COAD_genes_meta <- COAD_data_obj@rowRanges %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type, width) # Keep gene_id (Ensembl) as primary key

COAD_genes_meta$width_kb <- COAD_genes_meta$width / 1000

#####let's try GSEA to see if there are pathway differences ###########
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
    sample_ids %in% AK_samples_transcriptomic ~ "AK",
    sample_ids %in% AKP_samples_transcriptomic ~ "AKP",
    sample_ids %in% AKPS_samples_transcriptomic ~ "AKPS",
    TRUE ~ NA_character_
  )
) %>% filter(!is.na(group))

dds <- DESeqDataSetFromMatrix(countData = raw_counts_unique[, sample_metadata$sample_id],
                              colData = sample_metadata,
                              design = ~ group)

dds <- DESeq(dds)

comparisons <- list(
  c("AKP", "AK"),
  c("AKPS", "AKP"))

setwd("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Transcriptomic Analysis")

for (cmp in comparisons) {
  res <- results(dds, contrast = c("group", cmp[1], cmp[2]))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(COAD_genes_meta %>% dplyr::select(gene_id, gene_name), 
              by = c("ensembl_id" = "gene_id")) %>%
    filter(!is.na(stat) & !is.na(gene_name)) %>%
    arrange(desc(stat)) %>%
    dplyr::select(gene_name, stat)
  
  fname <- paste0(cmp[1], "_vs_", cmp[2], "_ranked_genes2.rnk")
  write.table(res_df, file = fname, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

AKvAKP_neg <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Transcriptomic Analysis/Output 2/AKPvAK2.GseaPreranked.1765669685143/gsea_report_for_na_neg_1765669685143.tsv")
AKvAKP_pos <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Transcriptomic Analysis/Output 2/AKPvAK2.GseaPreranked.1765669685143/gsea_report_for_na_pos_1765669685143.tsv")
AKvAKP <- rbind(AKvAKP_neg,AKvAKP_pos)
names(AKvAKP)[names(AKvAKP) == "NAME"] <- "pathway"
names(AKvAKP)[names(AKvAKP) == "FDR q-val"] <- "padj"
AKvAKP$NES <- as.numeric(AKvAKP$NES)

AKPvAKPS_neg <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Transcriptomic Analysis/Output 2/AKPSvAKP2.GseaPreranked.1765669769481/gsea_report_for_na_neg_1765669769481.tsv")
AKPvAKPS_pos <- read_tsv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/GSEA/Transcriptomic Analysis/Output 2/AKPSvAKP2.GseaPreranked.1765669769481/gsea_report_for_na_pos_1765669769481.tsv")
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
  
  write.xlsx(nes_matrix, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Final Data/Figure 8/8C_Heatmap/immmuneNES.xlsx")
  
  # 7. Plot full heatmap
  pheatmap(nes_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",
           display_numbers = plus_matrix,
           number_color = "black",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Hallmark Pathway Enrichment",
           fontsize_row = 20,
           fontsize_col = 20, 
           fontsize_number = 22,
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


AKvAKP <- AKvAKP[!is.na(AKvAKP$NES),]
AKPvAKPS <- AKPvAKPS[!is.na(AKPvAKPS$NES),]

gsea_results_list <- list(
  "AKP_vs_AK" = AKvAKP,
  "AKPS_vs_AKP" = AKPvAKPS
)

# Named vector to control column order and labels
comparisons <- c("AKP_vs_AK" = "AK vs. AKP",
                 "AKPS_vs_AKP"= "AKP vs. AKPS")

gsea_results_list_immmune_AKvAKP <- gsea_results_list$AKP_vs_AK[gsea_results_list$AKP_vs_AK$pathway %in% immune_pathways,]
gsea_results_list_oncogenic_AKvAKP <- gsea_results_list$AKP_vs_AK[gsea_results_list$AKP_vs_AK$pathway %in% oncogenic_pathways,]

gsea_results_list_immmune_AKPvAKPS <- gsea_results_list$AKPS_vs_AKP[gsea_results_list$AKPS_vs_AKP$pathway %in% immune_pathways,]
gsea_results_list_oncogenic_AKPvAKPS <- gsea_results_list$AKPS_vs_AKP[gsea_results_list$AKPS_vs_AKP$pathway %in% oncogenic_pathways,]



gsea_results_list_immune <- list(
  "AKP_vs_AK" = gsea_results_list_immmune_AKvAKP,
  "AKPS_vs_AKP" = gsea_results_list_immmune_AKPvAKPS
)


gsea_results_list_oncogenic <- list(
  "AKP_vs_AK" = gsea_results_list_oncogenic_AKvAKP,
  "AKPS_vs_AKP" = gsea_results_list_oncogenic_AKPvAKPS
)

########finally make the plots now ##########


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


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/hallMarkHeatMapHuman_immune_TRANSCRIPTOMIC_Dec13.pdf",hallMarkHeatMapHuman_IMMUNE, width = 18, height = 8, dpi=600)
ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/hallMarkHeatMapHuman_oncogenic_TRANSCRIPTOMIC_Dec13.pdf",hallMarkHeatMapHuman_ONCOGENIC, width = 18, height = 8, dpi=600)


###### now make volcano plots ######## 

# Function: shrinkage, adj p-values, top genes, and return ggplot object
make_volcano <- function(dds, contrast, plot_title, extra_labels = NULL,
                         fc_cut = 2, p_cut = 0.05, n_label = 10,
                         show_legend = TRUE) {
  
  # DESeq2 results with shrinkage
  res <- lfcShrink(dds, contrast = contrast, type = "ashr")
  
  # Add gene symbols and biotype (protein-coding filter)
  res <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(
      COAD_genes_meta %>% 
        dplyr::select(gene_id, gene_name, gene_type),
      by = c("ensembl_id" = "gene_id")
    ) %>%
    dplyr::filter(gene_type == "protein_coding")   # ✅ KEEP ONLY PROTEIN CODING
  
  # Replace rownames with HGNC symbols (required by EnhancedVolcano)
  rownames(res) <- make.unique(res$gene_name)
  
  # Adjust p-values (FDR)
  res$fdr_adj_pvalue <- p.adjust(res$pvalue, method = "fdr")
  
  # Helper to select gene labels (top sig up/down)
  get_top_genes <- function(df, fc_filter, n = 10) {
    subset_df <- df[fc_filter(df$log2FoldChange), ]
    subset_df <- subset_df[order(subset_df$fdr_adj_pvalue), ]
    subset_df$gene_name[1:min(n, nrow(subset_df))]
  }
  
  # Identify genes to label
  up_genes   <- get_top_genes(res, function(fc) fc >= fc_cut, n_label)
  down_genes <- get_top_genes(res, function(fc) fc <= -fc_cut, n_label)
  
  select_labels <- unique(c(up_genes, down_genes, extra_labels, "REG4"))
  
  
  # Volcano plot
  p <- EnhancedVolcano(
    res,
    lab = res$gene_name,
    selectLab = select_labels,
    x = "log2FoldChange",
    y = "fdr_adj_pvalue",
    title = plot_title,
    pCutoff = p_cut,
    FCcutoff = fc_cut,
    drawConnectors = TRUE,
    maxoverlapsConnectors = Inf,
    lengthConnectors = unit(0.0075, "npc"),
    boxedLabels = FALSE
  )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  
  return(p)
}

# Create plots
p1 <- make_volcano(dds, c("group", "AKP", "AK"), "AK vs AKP", show_legend = TRUE)
p2 <- make_volcano(dds, c("group", "AKPS", "AKP"), "AKP vs AKPS", show_legend = FALSE)

ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/volcanoAKvAKP_TRANSCRIPTOMIC_Dec13.pdf",p1, width = 10, height = 8)
ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/volcanoAKPvAKPS_TRANSCRIPTOMIC_Dec13.pdf",p2, width = 10, height = 8)


######## now make barplot illustrating immune-activation score for AK vs. AKP vs. AKPS ############# 

# Load required packages
library(dplyr)
library(ggplot2)
library(readr)

#------------------------------
# 1. Read in the CIBERSORT files
#------------------------------
ak <- read_csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Transcriptomic/CIBERSORT_AK_immunefraction.csv")
akp <- read_csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Transcriptomic/CIBERSORT_AKP_immunefraction.csv")
akps <- read_csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/output/Transcriptomic/CIBERSORT_AKPS_immunefraction.csv")

# Add group labels
ak$Group <- "AK"
akp$Group <- "AKP"
akps$Group <- "AKPS"

# Combine into one dataframe
df <- bind_rows(ak, akp, akps)

#------------------------------
# 2. Define immune populations to include
#------------------------------
immune_cells <- c(
  "T cells CD8",
  "T cells CD4 memory activated",
  "NK cells resting",
  "NK cells activated",
  "Macrophages M1",
  "Macrophages M2",
  "Plasma cells",
  "T cells regulatory (Tregs)",
  "Neutrophils"
)

# Check which of those are present in your file
immune_cells <- immune_cells[immune_cells %in% colnames(df)]
print(immune_cells)

#------------------------------
# 3. Compute immune activation score
#------------------------------
df <- df %>%
  mutate(immune_score = rowMeans(dplyr::select(., all_of(immune_cells)), na.rm = TRUE))

#------------------------------
# 4. Summarize per group with mean and 95% CI
#------------------------------
summary_df <- df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_score = mean(immune_score, na.rm = TRUE),
    sd = sd(immune_score, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_lower = mean_score - qt(0.975, df = n - 1) * se,
    ci_upper = mean_score + qt(0.975, df = n - 1) * se
  )

print(summary_df)

#------------------------------
# 5. Plot mean ± 95% CI barplot
#------------------------------
ggplot(summary_df, aes(x = Group, y = mean_score, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.8) +
  labs(
    title = "Immune Activation Score by Group",
    y = "Mean Immune Activation Score",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Optional: Save outputs
ggsave("immune_activation_barplot_Dec12.png", width = 5, height = 4, dpi = 300)
write_csv(summary_df, "immune_activation_stats.csv")

##### per Dr. Hwang's request create vendiagrams for volcano plots for mouse data ##### 

library(VennDiagram)
library(tidyverse)
library(ggVennDiagram)

# Define significance thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Convert all results to data frames with gene + stats
res_AKvAKP_df <- as.data.frame(results(dds, contrast = c("grp", "AKP_Tu", "AK_Tu")))
res_AKPvAKPS_df <- as.data.frame(results(dds, contrast = c("grp", "AKPS_Tu", "AKP_Tu")))
res_AKPSvLungMet_df <- as.data.frame(results(dds, contrast = c("grp", "AKPS_LungMet", "AKPS_Tu")))
res_AKPSvMet_df <- as.data.frame(results(dds, contrast = c("grp", "AKPS_Met", "AKPS_Tu")))
res_AKPSvTdLN_df <- as.data.frame(results(dds, contrast = c("grp", "AKPS_TdLN", "AKPS_Tu")))

# Add gene column
for (obj in c("res_AKvAKP_df", "res_AKPvAKPS_df", "res_AKPSvLungMet_df", "res_AKPSvMet_df", "res_AKPSvTdLN_df")) {
  tmp <- get(obj) %>% rownames_to_column("gene")
  assign(obj, tmp)
}

up_AK    <- res_AKvAKP_df       %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)
up_AKP   <- res_AKPvAKPS_df     %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)
up_AKPS  <- res_AKPSvMet_df     %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)
up_TdLN  <- res_AKPSvTdLN_df    %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)
up_Mets  <- res_AKPSvLungMet_df %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)

down_AK    <- res_AKvAKP_df       %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)
down_AKP   <- res_AKPvAKPS_df     %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)
down_AKPS  <- res_AKPSvMet_df     %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)
down_TdLN  <- res_AKPSvTdLN_df    %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)
down_Mets  <- res_AKPSvLungMet_df %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)

# --- combine lists ---
up_list <- list(AK = up_AK, AKP = up_AKP, AKPS = up_AKPS, TdLN = up_TdLN, Mets = up_Mets)
down_list <- list(AK = down_AK, AKP = down_AKP, AKPS = down_AKPS, TdLN = down_TdLN, Mets = down_Mets)

# --- 5-way Venn for UPREGULATED ---
venn_up <- venn.diagram(
  x = up_list,
  category.names = names(up_list),
  filename = "fiveWayVennDiagramUp.png",
  output = TRUE,
  fill = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3", "#377EB8"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = "black",
  main = "Upregulated Genes (padj < 0.05, LFC > 1)"
)
grid.newpage()
grid.draw(venn_up)

# --- 5-way Venn for DOWNREGULATED ---
venn_down <- venn.diagram(
  x = down_list,
  category.names = names(down_list),
  filename = "fiveWayVennDiagramDown.png",
  output = TRUE,
  fill = c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = "black",
  main = "Downregulated Genes (padj < 0.05, LFC < -1)"
)
grid.newpage()
grid.draw(venn_down)

### now look at CMS subtype progression across transcriptomic and genomic labels ########
### do genomic labels first ### 

CMS_classifications <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classifications_alluvial.csv")

CMS_classifications_barplot_grouped <- CMS_classifications[CMS_classifications$group %in% c("AK","AKP","AKPS"),]

CMS_classifications_barplot_grouped$RF <- replace(CMS_classifications_barplot_grouped$RF, is.na(CMS_classifications_barplot_grouped$RF), "Not Assigned")


CMS_classifications_barplot_grouped <- CMS_classifications_barplot_grouped %>%
  group_by(RF, group) %>%
  summarise(Freq =n(), .groups = "drop")

CMS_classifications_barplot_grouped <- CMS_classifications_barplot_grouped %>% group_by(group) %>% mutate(
  proportion = Freq/sum(Freq) 
) %>% ungroup()

#make alluvial based on genomic labeling

humanAlluvialCMS <- ggplot(CMS_classifications_barplot_grouped,
                           aes(x = group, stratum = RF, alluvium = RF,
                               y = proportion, fill = RF)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.7) +
  geom_stratum() +
  theme_minimal() +
  labs(title = "CMS Subtype Composition Across Tumor Progression",
       y = "", x = "Stage") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        legend.text  = element_text(size = 25),
        legend.title = element_text(size = 20)) + 
  guides(fill = guide_legend(title = NULL))


ggsave("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/Human_alluvialplot_CMS_Genomic_Dec10.pdf",humanAlluvialCMS, width = 14, height = 8)





#make alluvial based on transcriptomic labeling 


transcriptomic_samples <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CIBERSORT/Oct 12/lists/Transcriptomic/all_pred_results.csv")

transcriptomic_samples$Sample <- substr(transcriptomic_samples$Sample, 1, 15)

cms_labels <- read.csv("/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Human Data/CMS_classification_results.csv")

names(cms_labels)[names(cms_labels) == "X"] <- "Sample"

cms_labels$RF <- replace(cms_labels$RF, is.na(cms_labels$RF), "Not Assigned")

all_labeling <- merge(cms_labels ,transcriptomic_samples,by="Sample")

all_labeling <- all_labeling[,c("Sample","RF","Predicted_Class")]

names(all_labeling)[names(all_labeling) == "RF"] <- "CMS_subtype"

write.xlsx(all_labeling, "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/all_labeling.xlsx" )

####### now make alluvial plot 
df <- all_labeling %>%
  filter(Predicted_Class %in% c("AK", "AKP", "AKPS"))

# Count frequencies and compute proportions
df_prop_2 <- df %>%
  group_by(Predicted_Class, CMS_subtype) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(Predicted_Class) %>%
  mutate(proportion = Freq / sum(Freq)) %>%
  ungroup()

# Make alluvial plot
humanAlluvialCMS_transcriptomic <- ggplot(
  df_prop_2,
  aes(
    x = Predicted_Class,
    stratum = CMS_subtype,
    alluvium = CMS_subtype,
    y = proportion,
    fill = CMS_subtype
  )
) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.7) +
  geom_stratum() +
  theme_minimal() +
  labs(
    title = "CMS Subtype Composition Across Predicted Class",
    y = "",
    x = "Predicted Class"
  ) +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x  = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        legend.text  = element_text(size = 25),
        legend.title = element_text(size = 20)) +
  guides(fill = guide_legend(title = NULL))

# Save plot
ggsave(
  "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Output/Human_alluvialplot_CMS_Transcriptomic.pdf",
  humanAlluvialCMS_transcriptomic,
  width = 14,
  height = 8
)


