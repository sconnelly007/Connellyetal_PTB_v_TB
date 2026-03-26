library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidytext)
library(miRNAtap)
library(miRNAtap.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(tximeta)
library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(here)
library(DESeq2)
library(org.Hs.eg.db)
library(vroom)
library(RColorBrewer)
library(pheatmap)
library(UpSetR)
library(ComplexUpset)
library(patchwork)
library(EnhancedVolcano)
set.seed(42)

pca_plt <- function(dds_obj,
                    sample_prefix,nsub_s = 1000) {
  
  cols_cond <- c(Preterm = "#8C510A", Term = "#01665E")
  cols_race <- c(Caucasian = "#df9ed4",`African_American`="#469d76")
  #cols_age <- c("#ef8a47","#ffd06f","#aadce0","#72bcd5","#376795")
  #names(cols_age) <- levels(dds_obj$ga_bin)
  #cols_indication <- c(PTL   = "#a63603",  Term  = "#0f7ba2",  PPROM = "#fd8d3c")
  cols_sex <- c(Male="#374E55FF",Female= "#DF8F44FF")
  cols_bmi <- c("#f4a582","#4393c3","#762a83","#b2182b")
  names(cols_bmi) <- levels(dds_obj$bmi)
  
  vsd <- vst(dds_obj, blind = FALSE,nsub = nsub_s)
  
  pcaData <- plotPCA(vsd, intgroup = c("condition", "cb_match","race","ga_bin","indication_for_preterm_delivery","newborn_sex","bmi"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$SampleID <- rownames(pcaData)
  
  p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition,label = SampleID)) +
    geom_point(size =3) +
    ggrepel::geom_text_repel(size = 3,max.overlaps = 25) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    scale_color_manual(values=cols_cond)
  
  p2 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = race,label = SampleID)) +
    geom_point(size =3) +
    ggrepel::geom_text_repel(size = 3,max.overlaps = 25) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    scale_color_manual(values=cols_race)
  
  p3 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = newborn_sex,label = SampleID)) +
    geom_point(size =3) +
    ggrepel::geom_text_repel(size = 3,max.overlaps = 25) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    scale_color_manual(values=cols_sex)
  
  p4 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = bmi,label = SampleID)) +
    geom_point(size =3) +
    ggrepel::geom_text_repel(size = 3,max.overlaps = 25) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    scale_color_manual(values = cols_bmi)
  
  plt_m <- 
    p1 + p2 + p3 + p4 +
    plot_layout(ncol = 4,nrow = 1)
  ggsave(here(paste0("plots/",sample_prefix,"/",sample_prefix,"_pca.png")),plot = plt_m,width=20,height=20)
  return(list(vsd = vsd,
              plt = plt_m))
}

import_filter_mRNA <- function(txi_obj,
                               colD,
                               sample_prefix,
                               design_variables) {
  
  design_formula <- as.formula(
    paste("~", paste(design_variables, collapse = " + "))
  )
  
  dds <- DESeqDataSetFromTximport(
    txi = txi_obj,
    colData = colD,
    design = design_formula
  )
  
  #keep <- rowSums(counts(ddsMat)) > 0
  
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  
  # make a plot, distribution of gene counts
  gene_counts <- data.frame(gene_counts=rowSums(counts(dds)),keep=keep)
  
  # excludes 0 counts
  plot_filter <- ggplot(gene_counts, aes(x = log10(gene_counts + 1), fill = keep)) +
    geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("grey60", "steelblue"),
                      labels = c("Filtered", "Kept")) +
    labs(
      x = expression(log10("counts" + 1)),
      y = "Number of Genes",
      fill = ">=10 reads in ≥4 samples"
    ) +
    theme_minimal(base_size = 20)
  
  ggsave(here(paste0("plots/",sample_prefix,"/",sample_prefix,"_filtered_gene.png")),plot_filter,width=20,height=10,bg="white")
  
  # filter data
  dds <- dds[keep,]
  
  gene_ids <- rownames(dds)
  
  # remove Ensembl version suffix if present
  gene_ids_clean <- sub("\\..*$", "", gene_ids)
  
  rowData(dds)$ENSEMBL <- gene_ids_clean
  
  rowData(dds)$SYMBOL <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = gene_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first")
  
  return(list(dds = dds,
              plot_filter = plot_filter))
}

run_DE_analysis_mRNA <- function(txi_obj,
                                 colD,
                                 sample_prefix,
                                 design_variables,
                                 alpha_v_,
                                 log2fc_thresh_v,
                                 width=1500,height=2000) {
  print(paste0("Running: ",sample_prefix))
  # import + filter
  dds <- import_filter_mRNA(
    txi_obj = txi_obj,
    colD = colD,
    sample_prefix = sample_prefix,
    design_variables = design_variables
  )
  
  # sample distance heatmap
  sample_dist_plt <- samp_dist(dds$dds, sample_prefix)
  
  # PCA
  vsd <- pca_plt(dds$dds, sample_prefix)
  
  # differential expression
  de_results <- DE(dds$dds,
                   alpha_v = alpha_v_,
                   log2fc_thresh = log2fc_thresh_v)
  
  volcano <- volcano_plot(
    dds_obj = de_results$dds_de,
    res_obj  = de_results$res,
    sample_prefix = sample_prefix,
    alpha_v = alpha_v_,
    log2fc_thresh = log2fc_thresh_v,top_n = 2
  )
  
  de_plt <- plotDE_by_samp(
    results_df = de_results$res,
    vsd        = vsd$vsd,
    sample_prefix = sample_prefix,
    padj_thresh = alpha_v_,
    log2FC_t = log2fc_thresh_v,
    type = "mrna",width = width,height = height
  )
  
  return(list(dds = dds, 
              vsd = vsd, 
              de = de_results,
              samp_dist_plt = sample_dist_plt,
              volcano_plt = volcano,
              de_plt = de_plt
              ))
}

generate_upset_plt <- function(sigifiance_list,
                                    sample_name){
  
  custom_intersections <- list(
    "condition",
    c("condition","race"),
    c("condition","sex"),
    c("condition","bmi"),
    c("condition","race","bmi"),
    c("condition","sex","bmi"),
    c("condition","race","sex"),
    c("condition","race","sex","bmi")
  )
  
  incidence <- fromList(sigifiance_list)
  incidence_df <- as.data.frame(incidence)
  
  p <- ComplexUpset::upset(
    incidence_df,
    intersect = colnames(incidence_df),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        aes(fill = "Intersection")
      ) +
        scale_fill_manual(values = c("Intersection" = "black")) +
        ggtitle(paste0(sample_name," mRNA Analysis")) +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = "none"
        )
    ),
    set_sizes = (
      upset_set_size(
        aes(fill = "Set")
      ) +
        scale_fill_manual(values = c("Set" = "steelblue")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    )
  )
  
  # Save high-res PNG
  ggsave(
    filename = here(paste0("plots/compare_mRNA_upset_",sample_name,".png")),
    plot = p,
    width = 30,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  return(p)
}

# volcano
volcano_plot <- function(dds_obj,
                              res_obj,
                              sample_prefix,
                              alpha_v,
                              log2fc_thresh,
                              top_n,
                              type = "mrna"){
  plot_df <- as.data.frame(res_obj)
  rowdata_df <- as.data.frame(rowData(dds_obj))
  
  if (type == "mrna") {
    rowdata_df$common_ann <- ifelse(!is.na(rowData(dds_obj)$SYMBOL), rowData(dds_obj)$SYMBOL, rownames(dds_obj))
    rowdata_df$common_ann <- make.unique(rowdata_df$common_ann)
    rownames(plot_df) <- rowdata_df$common_ann
  }
  
  sig_genes <- plot_df %>% 
    filter(!is.na(padj) & padj < alpha_v & abs(log2FoldChange) > log2fc_thresh)
  
  sig_genes_up <- sig_genes %>%
    arrange(desc(log2FoldChange)) %>% 
    top_n(n = top_n) %>% 
    rownames()
  sig_genes_down <- sig_genes %>%
    arrange(log2FoldChange) %>% 
    top_n(n = top_n) %>% 
    rownames()
  genes_to_focus <- c(sig_genes_up,sig_genes_down)
  
  if (type == "mirna") {
    genes_to_focus <- rownames(plot_df)
  }
  
  plt <- EnhancedVolcano(plot_df,
                         lab = rownames(plot_df),
                         selectLab = genes_to_focus,
                         x = 'log2FoldChange',
                         y = 'padj',
                         xlab = bquote(~Log[2]~ 'fold change'),
                         ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                         title = "Differential Expression: mRNA",
                         subtitle = NULL,
                         caption = paste0("Total genes: ", nrow(plot_df)),
                         pCutoff = 0.05,
                         FCcutoff = 0.5,  
                         cutoffLineType = 'twodash',
                         cutoffLineWidth = 0.8,
                         pointSize = 2.5,
                         labSize = 4.5,
                         colAlpha = 0.9,
                         col = c("grey30", "forestgreen", "royalblue", "red2"),
                         legendLabels = c('Not sig.', 'Log2 FC', 'Adjusted p-value', 'P-value & Log2 FC'),
                         legendPosition = 'right',
                         legendLabSize = 12,
                         legendIconSize = 4.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.3,
                         max.overlaps = 15,
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE)
  png(filename = here(paste0("plots/",sample_prefix,"/",sample_prefix,"_volcano.png")),
      width=25,height=10,units="in",res=1200)
  print(plt)
  dev.off()
  return(plt)
}

get_sig_mRNA <- function(obj, alpha=0.05, lfc=0.5) {
  genes_id <- obj$de$res %>%
    as.data.frame() %>%
    filter(!is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc) %>%
    rownames()
  
  symbols <- rowData(obj$vsd) %>% 
    as.data.frame() %>%
    dplyr::filter(gene_id %in% genes_id) %>% 
    mutate(symbol_clean = ifelse(!is.na(SYMBOL), 
                                 SYMBOL, 
                                 gene_id)) %>% 
    pull(symbol_clean)
  return(list("ids" = genes_id,
              "symbols"=symbols))
}

plotDE_by_samp <- function(results_df,
                                vsd,
                                sample_prefix,
                                padj_thresh = 0.05,
                                log2FC_t = 0.5,
                                type = "mrna",
                           width=1500,height=2000) {
  
  cols_cond <- c(Preterm = "#8C510A", Term = "#01665E")
  cols_race <- c(Caucasian = "#df9ed4",`African_American`="#469d76")
  #cols_age <- c("#ef8a47","#ffd06f","#aadce0","#72bcd5","#376795")
  #names(cols_age) <- levels(vsd$ga_bin)
  #cols_indication <- c(PTL   = "#a63603",  Term  = "#0f7ba2",  PPROM = "#fd8d3c")
  cols_sex <- c(Male="#374E55FF",Female= "#DF8F44FF")
  cols_bmi <- c("#f4a582","#4393c3","#762a83","#b2182b")
  names(cols_bmi) <- levels(vsd$bmi)
  
  print(paste0("Are the rows the same? ",all(rownames(vsd) == rownames(results_df))))
  
  mat <- assay(vsd)
  if (type == "mrna"){
    rownames(mat) <- ifelse(!is.na(rowData(vsd)$SYMBOL), rowData(vsd)$SYMBOL, rownames(vsd))
    rownames(mat) <- sub("\\..*", "", rownames(mat))
    
    # filter out nont annotated genes
    is_annotated <- !str_detect(rownames(mat), "ENSG") 
    
    # DE genes
    keep <- !is.na(results_df$padj) &
      results_df$padj < padj_thresh &
      abs(results_df$log2FoldChange) > log2FC_t &
      is_annotated
  } else {
    keep <- !is.na(results_df$padj) &
      results_df$padj < padj_thresh &
      abs(results_df$log2FoldChange) > log2FC_t
  }
  
  mat <- mat[keep, , drop = FALSE]
  results_df_f <- results_df[keep, , drop = FALSE]
  
  # order by log2FC
  ord <- order(results_df_f$log2FoldChange,decreasing = TRUE)
  
  mat <- mat[ord, , drop = FALSE]
  
  # centering each gene's expression by subtracting mean expression of that gene
  mat <- mat - rowMeans(mat)
  
  anno_df <- data.frame(
    Condition = vsd$condition,
    #age = vsd$ga_bin,
    #race = vsd$race,
    row.names = colnames(vsd)
    #indication = vsd$indication_for_preterm_delivery,
    #sex = vsd$newborn_sex,
    #bmi = vsd$bmi
  )
  
  ann_colors = list(
    Condition = cols_cond
    #age=cols_age,
    #race=cols_race,
    #indication = cols_indication,
    #sex = cols_sex,
    #bmi = cols_bmi
  )
  
  conditions_ordered <- anno_df[order(anno_df$Condition, decreasing = TRUE), , drop = FALSE]
  mat <- mat[,rownames(conditions_ordered),drop = FALSE]
  
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#3B4CC0", "white", "#B40426")
  )
  
  ha <- HeatmapAnnotation(
    df = conditions_ordered,
    col = ann_colors,
    annotation_name_gp = gpar(fontface = "bold", fontsize = 20),
    annotation_legend_param = list(
      Condition = list(
        title_gp = gpar(fontface = "bold", fontsize = 24),
        labels_gp = gpar(fontsize = 20)
      )
    )
  )
  
  ht <- Heatmap(
    mat,
    name = "Centered log2FC",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 40),
    column_names_gp = gpar(fontsize = 40),
    row_names_side = "left",
    column_names_side = "bottom",
    row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 40)),
    column_names_max_height = unit(10, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.8),
    heatmap_legend_param = list(
      title = "Centered log2FC",
      title_gp = gpar(fontface = "bold", fontsize = 16),
      labels_gp = gpar(fontsize = 30),
      at = c(-2, -1, 0, 1, 2),
      legend_height = unit(4, "cm")
    )
  )
  png(here(paste0("plots/",sample_prefix,"/",sample_prefix,"_DE_genes_by_samp.png")),width=width,height=height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  return(ht)
}

# sample distance
samp_dist <- function(dds_obj,
                      sample_prefix,
                      nsub_s = 1000) {
  
  #sample distance
  #rld <- rlog(dds_obj, blind = FALSE)
  vsd <- vst(dds_obj, blind = FALSE,nsub=nsub_s)
  logcts <- assay(vsd)
  
  dists <- dist(t(logcts))
  dist_mat <- as.matrix(dists)
  
  colRamp <- viridis::viridis(100, option = "C")
  cols <- colRamp
  
  cols_cond <- c(Preterm = "#8C510A", Term = "#01665E")
  cols_race <- c(Caucasian = "#df9ed4",`African_American`="#469d76")
  #cols_age <- c("#ef8a47","#ffd06f","#aadce0","#72bcd5","#376795")
  #names(cols_age) <- levels(dds_obj$ga_bin)
  #cols_indication <- c(PTL   = "#a63603",  Term  = "#0f7ba2",  PPROM = "#fd8d3c")
  cols_sex <- c(Male="#374E55FF",Female= "#DF8F44FF")
  cols_bmi <- c("#f4a582","#4393c3","#762a83","#b2182b")
  names(cols_bmi) <- levels(dds_obj$bmi)
  
  anno_df <- data.frame(
    condition = dds_obj$condition,
    #age = dds_obj$ga_bin,
    race = dds_obj$race,
    row.names = colnames(dist_mat),
    #indication = dds_obj$indication_for_preterm_delivery,
    sex=dds_obj$newborn_sex,
    bmi=dds_obj$bmi
  )
  
  ann_colors = list(
    condition = cols_cond,
    #age=cols_age,
    race=cols_race,
    #indication = cols_indication,
    sex=cols_sex,
    bmi=cols_bmi
  )
  
  plot_ph <- pheatmap::pheatmap(dist_mat,
                                color = cols,
                                clustering_distance_rows = dists,
                                clustering_distance_cols = dists,
                                annotation_col = anno_df,
                                annotation_colors = ann_colors,
                                show_rownames = FALSE,
                                show_colnames = TRUE,
                                main=sample_prefix)
  
  save_pheatmap_png(plot_ph,
                    here(paste0("plots/",sample_prefix,"/",sample_prefix,"_sample_distance.png")),
                    width=2000,height=1500)
  return(plot_ph)
}

DE <- function(dds_obj,alpha_v,log2fc_thresh) {
  dds_de <- DESeq(dds_obj)
  
  res <- results(dds_de, alpha = alpha_v, lfcThreshold=log2fc_thresh, 
                    contrast = c("condition","Preterm","Term"),
                    altHypothesis="greaterAbs")
  
  resape <- DESeq2::lfcShrink(dds_de, coef = "condition_Preterm_vs_Term",
                              res =res,
                              type = "apeglm")
  return(list(dds_de = dds_de, res = res, res_shrunk = resape))
}

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

generate_upset_plt <- function(sigifiance_list,
                                    sample_name,
                                    type = "mRNA"){
  
  custom_intersections <- list(
    "condition",
    c("condition","race"),
    c("condition","sex"),
    c("condition","bmi"),
    c("condition","race","bmi"),
    c("condition","sex","bmi"),
    c("condition","race","sex"),
    c("condition","race","sex","bmi")
  )
  
  incidence <- fromList(sigifiance_list)
  incidence_df <- as.data.frame(incidence)
  
  p <- ComplexUpset::upset(
    incidence_df,
    intersect = colnames(incidence_df),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        aes(fill = "Intersection")
      ) +
        scale_fill_manual(values = c("Intersection" = "black")) +
        ggtitle(paste0(sample_name," ",type," Analysis")) +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = "none"
        )
    ),
    set_sizes = (
      upset_set_size(
        aes(fill = "Set")
      ) +
        scale_fill_manual(values = c("Set" = "steelblue")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    )
  )
  
  # Save high-res PNG
  ggsave(
    filename = here(paste0("plots/compare_",type,"_upset_",sample_name,".png")),
    plot = p,
    width = 30,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  return(p)
}


