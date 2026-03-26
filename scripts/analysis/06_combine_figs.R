library(cowplot)
library(here)
library(patchwork)
library(ComplexHeatmap)
library(tidyverse)
# load data ---------------------------------------------------------------
results_list_mRNA <- readRDS(here("data/processed_mrna/cb_results.rds"))

results_list_CB <- readRDS(here("data/processed_mirna/cb_mirna_results.rds"))

results_list_pl <- readRDS(here("data/processed_mirna/pl_results.rds"))

cb_upset <- readRDS(here("plots/cb_mirna_upset.rds"))
pl_upset <- readRDS(here("plots/pl_mirna_upset.rds"))
cb_mrna_upset <- readRDS(here("plots/cb_mrna_upset.rds"))

go_kegg_cb_up <- readRDS(here("data/processed_mirna/cb_go_res_up.rds"))
go_kegg_cb_down <- readRDS(here("data/processed_mirna/cb_go_res_down.rds"))

go_kegg_pl_up <- readRDS(here("data/processed_mirna/pl_go_res_up.rds"))
go_kegg_pl_down <- readRDS(here("data/processed_mirna/pl_go_res_down.rds"))

# Figure 1 ----------------------------------------------------------------
cb_mirna_de_plt <- results_list_CB$condition_race_sex_bmi$de_plt
cb_mirna_go_up <- go_kegg_cb_up$go_plt
cb_mirna_go_down <- go_kegg_cb_down$go_plt
heatmap_grob <- wrap_elements(
  grid.grabExpr(
    ComplexHeatmap::draw(
      cb_mirna_de_plt, 
      heatmap_legend_side = "right", 
      annotation_legend_side = "right",
      merge_legend = TRUE,
      padding = unit(c(1, 1, 1, 3), "cm")
    )
  )
)
combined_figure <- heatmap_grob + (cb_mirna_go_up / cb_mirna_go_down) +
  plot_layout(widths = c(2.5, 1)) +
  plot_annotation(tag_levels = 'A')
ragg::agg_png(
  filename = here("plots/final_figs/fig1.png"),
  width = 22,
  height = 10,
  units = "in",
  res = 600)
print(combined_figure)
dev.off()

# save indiv for poster 
ragg::agg_png(
  filename = here("plots/final_figs/fig1A.png"),
  width = 18,
  height = 20,
  units = "in",
  res = 600)
print(heatmap_grob)
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig1B.png"),
  width = 16,
  height = 12,
  units = "in",
  res = 600)
print(cb_mirna_go_up)
dev.off()


ragg::agg_png(
  filename = here("plots/final_figs/fig1C.png"),
  width = 16,
  height = 12,
  units = "in",
  res = 600)
print(cb_mirna_go_down)
dev.off()

# Figure 2 ----------------------------------------------------------------
pl_mirna_de_plt <- results_list_pl$condition_race_sex_bmi$de_plt
pl_mirna_go_up <- go_kegg_pl_up$go_plt
pl_mirna_go_down <- go_kegg_pl_down$go_plt
heatmap_grob <- wrap_elements(
  grid.grabExpr(
    ComplexHeatmap::draw(
      pl_mirna_de_plt, 
      heatmap_legend_side = "right", 
      annotation_legend_side = "right",
      merge_legend = TRUE,
      padding = unit(c(1, 1, 1, 3), "cm")
    )
  )
)
combined_figure <- heatmap_grob + (pl_mirna_go_up / pl_mirna_go_down) +
  plot_layout(widths = c(2.5, 1)) +
  plot_annotation(tag_levels = 'A')
ragg::agg_png(
  filename = here("plots/final_figs/fig2.png"),
  width = 22,
  height = 10,
  units = "in",
  res = 600)
print(combined_figure)
dev.off()

# save indiv for poster 
ragg::agg_png(
  filename = here("plots/final_figs/fig2A.png"),
  width = 18,
  height = 20,
  units = "in",
  res = 600)
print(heatmap_grob)
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig2B.png"),
  width = 16,
  height = 12,
  units = "in",
  res = 600)
print(pl_mirna_go_up)
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig2C.png"),
  width = 16,
  height = 12,
  units = "in",
  res = 600)
print(pl_mirna_go_down)
dev.off()

# Figure 3 ----------------------------------------------------------------
mrna_de_plt <- results_list_mRNA$condition_race_sex_bmi$de_plt
ht_up <- readRDS(here("results/ht_compare_miRNA_up_mRNA_down.rds"))
ht_down <- readRDS(here("results/ht_compare_miRNA_down_mRNA_up.rds"))

heatmap_grob <- wrap_elements(
  grid.grabExpr(
    ComplexHeatmap::draw(
      mrna_de_plt, 
      heatmap_legend_side = "right", 
      annotation_legend_side = "right",
      merge_legend = TRUE,
      padding = unit(c(1, 1, 1, 3), "cm")
    )
  )
)

heatmap_grob_ht_up <- wrap_elements(
  grid.grabExpr(
    ComplexHeatmap::draw(
      ht_up
    )
  )
)

heatmap_grob_ht_down <- wrap_elements(
  grid.grabExpr(
    ComplexHeatmap::draw(
      ht_down
    )
  )
)


ragg::agg_png(
  filename = here("plots/final_figs/fig3_A.png"),
  width = 24,
  height = 32,
  units = "in",
  res = 600)
heatmap_grob +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(size = 36, face = "bold")))
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig3_B.png"),
  width = 24,
  height = 10,
  units = "in",
  res = 600)
(heatmap_grob_ht_up + heatmap_grob_ht_down)+
  plot_annotation(tag_levels = list(c('B', 'C')),
                  theme = theme(plot.tag = element_text(size = 36, face = "bold")))
dev.off()

# indiv plotting
ragg::agg_png(
  filename = here("plots/final_figs/fig3_A_poster.png"),
  width = 24,
  height = 34,
  units = "in",
  res = 600)
print(heatmap_grob)
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig3_B_poster.png"),
  width = 30,
  height = 10,
  units = "in",
  res = 600)
print(heatmap_grob_ht_up)
dev.off()

ragg::agg_png(
  filename = here("plots/final_figs/fig3_C_poster.png"),
  width = 20,
  height = 20,
  units = "in",
  res = 600)
print(heatmap_grob_ht_down)
dev.off()

# supplementary figure 1 --------------------------------------------------
cb_mirna_gf <- results_list_CB$condition_race_sex_bmi$dds$plot_filter

pl_mirna_gf <- results_list_pl$condition_race_sex_bmi$dds$plot_filter

cb_mrna_gf <- results_list_mRNA$condition_race_sex_bmi$dds$plot_filter

ragg::agg_png(
  filename = here("plots/final_figs/combined_genefilter_highres.png"),
  width = 20,
  height = 20,
  units = "in",
  res = 600
)
cb_mirna_gf / pl_mirna_gf / cb_mrna_gf + plot_annotation(tag_levels = c("A","B","C"))
dev.off()

# supplementary figure 2 --------------------------------------------------
cb_mirna_pca <- results_list_CB$condition_race_sex_bmi$vsd$plt

pl_mirna_pca <- results_list_pl$condition_race_sex_bmi$vsd$plt

cb_mrna_pca <- results_list_mRNA$condition_race_sex_bmi$vsd$plt

combined_pca <- 
  wrap_elements(cb_mirna_pca) /
  wrap_elements(pl_mirna_pca) /
  wrap_elements(cb_mrna_pca) +
  plot_annotation(tag_levels = "A")

ragg::agg_png(
  filename = here("plots/final_figs/combined_PCA_highres.png"),
  width = 30,
  height = 20,
  units = "in",
  res = 600
)
combined_pca
dev.off()

# supplementary figure 3 --------------------------------------------------
combined_upset <- plot_grid(
  cb_upset,
  pl_upset,
  cb_mrna_upset,
  labels = c("A", "B","C"),
  ncol = 1,
  nrow = 3,
  label_size = 18
)

ragg::agg_png(
  filename = here("plots/final_figs/combined_upset_highres.png"),
  width = 20,
  height = 20,
  units = "in",
  res = 600
)
print(combined_upset)
dev.off()
