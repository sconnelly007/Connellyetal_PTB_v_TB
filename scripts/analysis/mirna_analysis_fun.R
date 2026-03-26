import_filter_mirna <- function(count_mat, 
                                coldata,
                                sample_prefix,
                                design_variables) {
  
  design_formula <- as.formula(
    paste("~", paste(design_variables, collapse = " + "))
  )
  
  ddsMat <- DESeqDataSetFromMatrix(countData = count_mat,
                                   colData = coldata,
                                   design = design_formula)

  smallestGroupSize <- 4
  keep <- rowSums(counts(ddsMat) >= 10) >= smallestGroupSize
  
  # make a plot, distribution of gene counts
  gene_counts <- data.frame(gene_counts=rowSums(counts(ddsMat)),keep=keep)
  
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
  ddsMat <- ddsMat[keep,]
  
  return(list(dds = ddsMat,
              plot_filter = plot_filter))
}


get_sig <- function(res, alpha, lfc) {
  res %>%
    as.data.frame() %>%
    filter(!is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc) %>%
    rownames()
}

run_DE_analysis_mirna <- function(counts, 
                                  coldata, 
                                  sample_prefix, 
                                  design_variables,
                                  alpha_v_,
                                  log2fc_thresh_v,
                                  nsub_s = 100) {
  print(paste0("Running: ",sample_prefix))
  # import + filter
  dds <- import_filter_mirna(
    count_mat = counts,
    coldata   = coldata,
    sample_prefix = sample_prefix,
    design_variables = design_variables
  )
  
  # sample distance heatmap
  sample_dist_plt <- samp_dist(dds$dds, sample_prefix,nsub_s = nsub_s)
  
  # PCA
  vsd <- pca_plt(dds$dds, sample_prefix,nsub_s = nsub_s)
  
  # differential expression
  de_results <- DE(dds$dds,alpha_v = alpha_v_,log2fc_thresh = log2fc_thresh_v)
  
  sig_genes <- get_sig(de_results$res,alpha = alpha_v_,lfc = log2fc_thresh_v)
  
  if (length(sig_genes) < 1) {
    message(sample_prefix, ": No significant DE genes found. Skipping downstream plots.")
    
    return(list(
      dds = dds,
      vsd = vsd,
      de  = de_results,
      sig = character(0)
    ))
  } else {
    volcano <- volcano_plot(
      dds_obj = de_results$dds_de,
      res_obj  = de_results$res,
      sample_prefix = sample_prefix,
      alpha_v = alpha_v_,
      log2fc_thresh = log2fc_thresh_v,
      top_n = 2,
      type = "mirna"
    )
    
    de_plt <- plotDE_by_samp(
      results_df = de_results$res,
      vsd        = vsd$vsd,
      sample_prefix = sample_prefix,
      padj_thresh = alpha_v_,
      log2FC_t = log2fc_thresh_v,
      type = "mirna",width = 1000,height = 1000
    )
    
    return(list(dds = dds, 
                vsd = vsd, 
                de = de_results,
                samp_dist_plt = sample_dist_plt,
                volcano_plt = volcano,
                de_plt = de_plt,
                sig = sig_genes
    ))
  }
}
