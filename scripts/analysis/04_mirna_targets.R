library(miRNAtap)
library(miRNAtap.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)

# miRNAtap to find miRNA targets in CB ------------------------------------------
results_list_CB <- readRDS(here("data/processed_mirna/cb_mirna_results.rds"))
alpha <- 0.05
log2FC_thresh <- 0.3
sig_up_cb <- results_list_CB$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange > log2FC_thresh) %>% 
  rownames()
sig_down_cb <- results_list_CB$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange < -log2FC_thresh) %>% 
  rownames()

# GO Enrichment in CB --------------------------------------------------
plot_miRNA_GO <- function(go_res_df,color) {
  # clean df
  plot_data <- go_res_df %>%
    mutate(
      # reformat p values
      KS_clean = as.numeric(gsub("< ", "", KS)), 
      # negative log of the p values
      Log10P = -log10(KS_clean),
      # wrap GO terms
      Term = str_wrap(Term, width = 45), 
      # reorder terms by significance
      Term = reorder(Term, Log10P)
    )
  ggplot(plot_data, aes(x = Log10P, y = Term)) +
    geom_col(show.legend = FALSE, width = 0.7,fill = color) +
    scale_y_reordered() +
    theme_minimal(base_size = 40) +
    # Labels
    labs(
      x = expression("-log"[10]*"(p-value)"),
      y = NULL
    ) #+
    # Fine-tuning the look
    # theme(
    #   strip.text = element_text(face = "bold", size = 12),
    #   panel.grid.major.y = element_blank(),
    #   plot.title = element_text(face = "bold"),
    #   axis.text.y = element_text(size = 10)
    # )
}
miRNA_target_GO <- function(significant_miRNAs,
                            ontology = c("BP","CC","MF"),
                            top_n,
                            go_thresh = 0.05,
                            kegg_p_cutoff = 0.05,
                            col = "#e76254"){
  set.seed(1234)
  all_preds <- map(significant_miRNAs, function(m) {
    pred <- getPredictedTargets(m, species = "hsa", method = "geom", min_src = 2)
    
    if (is.null(pred) || nrow(pred) == 0) return(NULL)
    
    pred %>%
      as.data.frame() %>%
      rownames_to_column("entrez") %>%
      mutate(miRNA = m)
  }) %>%
    bind_rows()
  
  combined_ranked_cb <- all_preds %>%
    group_by(entrez) %>%
    summarise(rank_product = min(rank_product, na.rm = TRUE), .groups = "drop")
  
  rankedGenes <- combined_ranked_cb$rank_product
  names(rankedGenes) <- combined_ranked_cb$entrez
  
  selection <- function(x) TRUE
  
  all_ont <- map(ontology,function(ont){
    allGO2genes <- annFUN.org(whichOnto=ont, feasibleGenes = NULL,mapping="org.Hs.eg.db", ID = "entrez")
    GOdata <- new('topGOdata', ontology = ont, 
                  allGenes = rankedGenes,
                  annot = annFUN.GO2genes, 
                  GO2genes = allGO2genes,
                  geneSel = selection, 
                  nodeSize=10)
    results.ks <- runTest(GOdata, 
                          algorithm = "classic", 
                          statistic = "ks")
    allRes <- GenTable(GOdata, 
                       KS = results.ks, 
                       orderBy = "KS",
                       topNodes = length(usedGO(GOdata)))
    
    # select top n for each ontology
    allRes %>% 
      dplyr::mutate(
        KS = as.numeric(KS),
        Term = sapply(`GO.ID`, function(go) Term(GOTERM[[go]])),
        ont = ont) %>%
      filter(KS < go_thresh) %>%
      dplyr::arrange(KS)%>%
      slice_head(n = top_n)
    }) %>%
    bind_rows()
  
  # adjust KEGG mapping so that the genes are ranked in descending order
  # rank_kegg <- combined_ranked_cb %>% 
  #   arrange(desc(rank_product))
  # 
  # # pull vector
  # ranked_list <- rank_kegg$rank_product
  # names(ranked_list) <- rank_kegg$entrez
  # 
  # gse_KEGG_res <- gseKEGG(
  #   geneList = ranked_list,
  #   nPermSimple = 10000,
  #   minGSSize    = 3,
  #   maxGSSize    = 800,
  #   pvalueCutoff = kegg_p_cutoff,
  #   pAdjustMethod = "BH")
  
  go_plt <- plot_miRNA_GO(all_ont,col)
  
  return(list(targets = all_preds,
         rank_v = rankedGenes,
         GO_res = all_ont,
         go_plt = go_plt#,
         #kegg_p = gse_KEGG_res
         ))
}


go_kegg_cb_up <- miRNA_target_GO(sig_up_cb,
                                 ontology = c("BP"),
                                 top_n = 10,
                                 go_thresh = 0.05,col = "#e76254")
                                 
saveRDS(go_kegg_cb_up,here("data/processed_mirna/cb_go_res_up.rds"))
#go_kegg_cb_up$go_plt
go_kegg_cb_down <- miRNA_target_GO(sig_down_cb,
                                   ontology = c("BP"),
                                   top_n = 10,
                                   go_thresh = 0.05,col = "#e76254")
#go_kegg_cb_down$go_plt
saveRDS(go_kegg_cb_down,here("data/processed_mirna/cb_go_res_down.rds"))

# miRNAtap to find miRNA targets in PL ------------------------------------
results_list_pl <- readRDS(here("data/processed_mirna/pl_results.rds"))
alpha <- 0.05
log2FC_thresh <- 0.3
sig_up_pl <- results_list_pl$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange > log2FC_thresh) %>% 
  rownames()
sig_down_pl <- results_list_pl$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange < -log2FC_thresh) %>% 
  rownames()

go_kegg_pl_up <- miRNA_target_GO(sig_up_pl,
                                 ontology = c("BP"),
                                 top_n = 10,
                                 go_thresh = 0.05,col = "#72bcd5")
saveRDS(go_kegg_pl_up,here("data/processed_mirna/pl_go_res_up.rds"))
go_kegg_pl_up$go_plt

go_kegg_pl_down <- miRNA_target_GO(sig_down_pl,
                                   ontology = c("BP"),
                                   top_n = 10,
                                   go_thresh = 0.05,col = "#72bcd5")
saveRDS(go_kegg_pl_down,here("data/processed_mirna/pl_go_res_down.rds"))
go_kegg_pl_down$go_plt

# do enrichment in significance mRNA --------------------------------------
results_list_mRNA <- readRDS(here("data/processed_mrna/cb_results.rds"))
alpha <- 0.05
log2FC_thresh <- 0.3

results_list_mRNA_res <- results_list_mRNA$condition_race_sex_bmi$de$res
vsd <- results_list_mRNA$condition_race_sex_bmi$vsd$vsd
results_list_mRNA_res$SYMBOL <- ifelse(!is.na(rowData(vsd)$SYMBOL), rowData(vsd)$SYMBOL, rownames(vsd))
results_list_mRNA_res$gene_ids <- sub("\\..*", "", rownames(results_list_mRNA_res))

sig_up_cb_mrna <- results_list_mRNA_res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange > log2FC_thresh) %>% 
  pull(gene_ids)
sig_down_cb_mrna <- results_list_mRNA_res %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < alpha & log2FoldChange < -log2FC_thresh) %>% 
  pull(gene_ids)
universe_genes <- results_list_mRNA_res %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(padj)) %>% 
  pull(gene_ids)

# test genes that are upregulated
GO_results_cb_mrna_up <- enrichGO(gene = sig_up_cb_mrna, 
                                  universe      = universe_genes,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = "ENSEMBL", 
                                  ont = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = alpha,
                                  qvalueCutoff = 0.2,
                                  readable = TRUE) 

GO_results_cb_mrna_up_df <- GO_results_cb_mrna_up %>% 
  as.data.frame() %>% 
  group_by(geneID) %>% 
  slice_head(n = 2)

# test genes that are downregulated
GO_results_cb_mrna_down <- enrichGO(gene = sig_down_cb_mrna, 
                                    universe      = universe_genes,
                                    OrgDb = "org.Hs.eg.db", 
                                    keyType = "ENSEMBL", 
                                    ont = "BP",
                                    pvalueCutoff = alpha,
                                    qvalueCutoff = 0.2,
                                    readable = TRUE)
GO_results_cb_mrna_down_df <- GO_results_cb_mrna_down %>% 
  as.data.frame() %>% 
  group_by(geneID) %>% 
  slice_head(n = 2)

# return targets of all up or down mirnas ---------------------------------
up_cb <- results_list_CB$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(log2FoldChange > log2FC_thresh) %>% 
  rownames()

all_preds_up <- map(up_cb, function(m) {
  pred <- getPredictedTargets(m, species = "hsa", method = "geom", min_src = 2)
  
  if (is.null(pred) || nrow(pred) == 0) return(NULL)
  
  pred %>%
    as.data.frame() %>%
    rownames_to_column("entrez") %>%
    mutate(miRNA = m)
}) %>%
  bind_rows()
saveRDS(all_preds_up,here("data/processed_mirna/cb_targets_up.rds"))

down_cb <- results_list_CB$condition_race_sex_bmi$de$res %>% 
  as.data.frame() %>% 
  dplyr::filter(log2FoldChange < -log2FC_thresh) %>% 
  rownames()

all_preds_down <- map(down_cb, function(m) {
  pred <- getPredictedTargets(m, species = "hsa", method = "geom", min_src = 2)
  
  if (is.null(pred) || nrow(pred) == 0) return(NULL)
  
  pred %>%
    as.data.frame() %>%
    rownames_to_column("entrez") %>%
    mutate(miRNA = m)
}) %>%
  bind_rows()
saveRDS(all_preds_down,here("data/processed_mirna/cb_targets_down.rds"))

