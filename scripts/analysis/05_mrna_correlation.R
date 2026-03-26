alpha <- 0.05
log2fc_thresh <- 0.3
library(here)
library(miRNAtap)
library(miRNAtap.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

##### Correlation 
cor_one <- function(mi, mi_df,gene,mrna_df) {
  x <- as.numeric(mi_df[mi, ])
  y <- as.numeric(mrna_df[gene, ])
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  tibble(miRNA = mi, target_id = gene, r = unname(ct$estimate), pval = ct$p.value, numb = length(x))
}
cor_pairs <- function(pairs_df, 
                      mi_df,
                      mrna_df,
                      direction_label) {
  map(
    seq_len(nrow(pairs_df)),
    function(i) {
      mi   <- pairs_df[["miRNA"]][i]
      gene <- pairs_df[["target_ensembl"]][i]
      sym  <- pairs_df[["SYMBOL"]][i]
      
      out <- cor_one(mi, mi_df, gene, mrna_df)
      out$target_symbol <- sym
      out$direction     <- direction_label
      out
    }
  ) %>%
    list_rbind() %>% 
    mutate(padj = p.adjust(pval, method = "BH"))
}

plot_interaction_heatmap <- function(cor_df,
                                     file_out,
                                     title,
                                     bg = "white",
                                     width,
                                     height) {
  
  # wide matrix of Pearson r (rows = miRNA, cols = mRNA)
  mat <- cor_df %>%
    select(miRNA, target_symbol, r) %>%
    pivot_wider(names_from = target_symbol, values_from = r, values_fill = NA) %>%
    tibble::column_to_rownames("miRNA") %>%
    as.matrix()
  
  # number of DE targets interactions per miRNA
  row_deg <- rowSums(!is.na(mat))
  row_max <- max(row_deg, na.rm = TRUE)
  row_breaks <- unique(round(seq(0, row_max, length.out = 3)))
  
  # number of DE miRNA interactions per gene
  col_deg <- colSums(!is.na(mat))
  col_max <- max(col_deg, na.rm = TRUE)
  col_breaks <- unique(round(seq(0, col_max, length.out = 3)))
  
  mat[is.na(mat)] <- 0
  
  row_ha <- rowAnnotation(
    `Targets` = anno_barplot(row_deg, 
                             gp = gpar(fill = "#72bcd5",col = NA,fontsize = 40),
                             axis_param = list(at = row_breaks, labels = row_breaks),
                             bar_width = 0.7),
    annotation_name_side = "bottom"
  )
  top_ha <- HeatmapAnnotation(
    `miRNAs` = anno_barplot(col_deg, 
                            gp = gpar(fill = "#ef8a47", 
                                      col = NA,fontsize = 40), 
                            axis_param = list(at = col_breaks, labels = col_breaks),
                            bar_width = 0.7),
    annotation_name_side = "left"
  )
  
  col_fun <- colorRamp2(c(-1, -0.5), c("#2166ac", "#f0f0f0"))
  
  ht <- Heatmap(
    mat,
    name = "Pearson r",
    col = col_fun,
    na_col = "#f0f0f0",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE, 
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 40),
    column_names_gp = gpar(fontsize = 40),
    row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 40)),
    column_names_max_height = unit(10, "cm"),
    row_names_side = "left", 
    column_names_side = "bottom",
    right_annotation = row_ha,
    top_annotation = top_ha,
    heatmap_legend_param = list(title = "Pearson r", at = c(-1,-0.5))
  )
  
  # save
  png(file_out <- file_out, width = width, height = height, res = 300, bg = bg)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title=title,
       column_title_gp=grid::gpar(fontsize=16))
  dev.off()
  return(ht)
}


# load data ---------------------------------------------------------------
results_list_mRNA <- readRDS(here("data/processed_mrna/cb_results.rds"))
results_list_CB <- readRDS(here("data/processed_mirna/cb_mirna_results.rds"))
go_kegg_cb_up <- readRDS(here("data/processed_mirna/cb_go_res_up.rds"))
go_kegg_cb_down <- readRDS(here("data/processed_mirna/cb_go_res_down.rds"))

mrna_vsd_mat <- assay(results_list_mRNA$condition_race_sex_bmi$vsd$vsd)
rownames(mrna_vsd_mat) <- sub("\\..*", "",rownames(mrna_vsd_mat))
mirna_vsd_mat <- assay(results_list_CB$condition_race_sex_bmi$vsd$vsd) 
colnames(mirna_vsd_mat) <- gsub("\\.", "-", colnames(mirna_vsd_mat))

common_samples <- intersect(colnames(mirna_vsd_mat), colnames(mrna_vsd_mat))
mirna_vsd_mat <- mirna_vsd_mat[, common_samples, drop = FALSE]
mrna_vsd_mat <- mrna_vsd_mat[, common_samples, drop = FALSE]

symb_table <- data.frame(gene_id = rownames(mrna_vsd_mat))
symb_table$SYMBOL <- ifelse(!is.na(rowData(results_list_mRNA$condition_race_sex_bmi$vsd$vsd)$SYMBOL), 
                            rowData(results_list_mRNA$condition_race_sex_bmi$vsd$vsd)$SYMBOL, 
                            rownames(mrna_vsd_mat))

# compare all up miRNA and down mRNA --------------------------------------
genes_pass_down_mrna <- results_list_mRNA$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  filter(log2FoldChange < -log2fc_thresh) %>% 
  rownames_to_column("gene") %>% 
  mutate(gene = sub("\\..*", "",gene)) %>% 
  pull(gene)
mrna_vsd_mat_down <- mrna_vsd_mat[genes_pass_down_mrna,,drop= F]

all_preds_up <- readRDS(here("data/processed_mirna/cb_targets_up.rds"))
ids <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(all_preds_up$entrez), 
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENTREZID"
) %>% 
  distinct()

master_targets <- all_preds_up %>% 
  mutate(entrez = as.character(entrez)) %>% 
  left_join(ids, by = c("entrez" = "ENTREZID")) %>% 
  filter(ENSEMBL %in% rownames(mrna_vsd_mat_down)) %>%
  rename(target_ensembl = ENSEMBL)

mirnatap_correlations_all_up <- cor_pairs(
  pairs_df = master_targets,
  mi_df = mirna_vsd_mat,
  mrna_df = mrna_vsd_mat_down,
  direction_label = "miRNA_up_mRNA_down"
)
saveRDS(mirnatap_correlations_all_up,here("results/miRNA_up_mRNA_down_all.rds"))

mirnatap_correlations_all_up_sig <- mirnatap_correlations_all_up %>% 
  filter(padj < 0.01,r <= -0.5) %>% 
  mutate(id_symbol = paste0(target_symbol,"_",target_id)) %>% 
  add_count(miRNA) %>% 
  filter(n > 2) %>% 
  select(-n)

ht_up <- plot_interaction_heatmap(mirnatap_correlations_all_up_sig,
                         file_out = here("plots/compare_miRNA_up_mRNA_down.png"),
                         title = "miRNA upreg. - mRNA downreg.",
                         width=3500,
                         height=3000)
base::saveRDS(ht_up,here("results/ht_compare_miRNA_up_mRNA_down.rds"))


# GO Analysis!
genes_blk1_up <- mirnatap_correlations_all_up_sig %>% 
  filter(miRNA == "hsa-miR-2110") %>% 
  pull(target_id)

GO_results_mirna_up_mRNA_down_blk1 <- enrichGO(gene = genes_blk1_up, 
                                               universe      = rownames(mrna_vsd_mat_down),
                                               OrgDb = "org.Hs.eg.db", 
                                               keyType = "ENSEMBL", 
                                               ont = "BP",
                                               pvalueCutoff = 1,
                                               qvalueCutoff = 1,
                                               readable = TRUE)



# compare all up mrna and down mirna --------------------------------------
genes_pass_up_mrna <- results_list_mRNA$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  filter(log2FoldChange > log2fc_thresh) %>% 
  rownames_to_column("gene") %>% 
  mutate(gene = sub("\\..*", "",gene)) %>% 
  pull(gene)

mrna_vsd_mat_up <- mrna_vsd_mat[genes_pass_up_mrna,,drop= F]

all_preds_down <- readRDS(here("data/processed_mirna/cb_targets_down.rds"))

ids <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(all_preds_down$entrez), 
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENTREZID"
) %>% 
  distinct()

master_targets <- all_preds_down %>% 
  mutate(entrez = as.character(entrez)) %>% 
  left_join(ids, by = c("entrez" = "ENTREZID")) %>% 
  filter(ENSEMBL %in% rownames(mrna_vsd_mat_up)) %>%
  rename(target_ensembl = ENSEMBL)

mirnatap_correlations_down <- cor_pairs(
  pairs_df = master_targets,
  mi_df = mirna_vsd_mat,
  mrna_df = mrna_vsd_mat_up,
  direction_label = "miRNA_down_mRNA_up"
)
saveRDS(mirnatap_correlations_down,here("results/miRNA_up_mRNA_down_all.rds"))

mirnatap_correlations_all_down_sig <- mirnatap_correlations_down %>% 
  filter(padj < 0.01,r <= -0.5) %>% 
  mutate(id_symbol = paste0(target_symbol,"_",target_id)) %>% 
  add_count(miRNA) %>% 
  filter(n > 2) %>% 
  select(-n)

ht_down <- plot_interaction_heatmap(mirnatap_correlations_all_down_sig,
                         file_out = here("plots/compare_miRNA_down_mRNA_up.png"),
                         title = "miRNA downreg. - mRNA upreg.",
                         width=3500,
                         height=3000)
base::saveRDS(ht_down,here("results/ht_compare_miRNA_down_mRNA_up.rds"))

# GO Analysis!

# do down regulated mrna and upregulated sig mirna --------------------------------------------
# up_targets <- go_kegg_cb_up$targets
# 
# ids <- AnnotationDbi::select(
#   org.Hs.eg.db,
#   keys = as.character(up_targets$entrez), 
#   columns = c("ENSEMBL", "SYMBOL"),
#   keytype = "ENTREZID"
# ) %>% 
#   distinct()
# 
# master_targets <- up_targets %>% 
#   mutate(entrez = as.character(entrez)) %>% 
#   left_join(ids, by = c("entrez" = "ENTREZID")) %>% 
#   filter(ENSEMBL %in% rownames(mrna_vsd_mat_down)) %>%
#   rename(target_ensembl = ENSEMBL)
# 
# mirnatap_correlations <- cor_pairs(
#   pairs_df = master_targets,
#   mi_df = mirna_vsd_mat,
#   mrna_df = mrna_vsd_mat_down,
#   direction_label = "miRNA_up_mRNA_down"
# )
# 
# 
# # correlate with up mrna and down sig mirna -----------------------------------
# down_targets <- go_kegg_cb_down$targets
# 
# ids <- AnnotationDbi::select(
#   org.Hs.eg.db,
#   keys = as.character(down_targets$entrez), 
#   columns = c("ENSEMBL", "SYMBOL"),
#   keytype = "ENTREZID"
# ) %>% 
#   distinct()
# 
# master_targets <- down_targets %>% 
#   mutate(entrez = as.character(entrez)) %>% 
#   left_join(ids, by = c("entrez" = "ENTREZID")) %>% 
#   filter(ENSEMBL %in% rownames(mrna_vsd_mat_up)) %>%
#   rename(target_ensembl = ENSEMBL)
# 
# # split dataframe by mirna
# mirna_split_list <- split(down_targets, down_targets$miRNA)
# 
# mirnatap_correlations <- cor_pairs(
#   pairs_df = master_targets,
#   mi_df = mirna_vsd_mat,
#   mrna_df = mrna_vsd_mat_up,
#   direction_label = "miRNA_down_mRNA_up"
# )
# 
