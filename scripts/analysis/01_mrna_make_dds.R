library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidytext)

designs <- list(
  condition        = c("condition"),
  condition_race   = c("condition","race"),
  condition_sex    = c("condition","newborn_sex"),
  condition_bmi   = c("condition","bmi_raw"),
  condition_race_bmi   = c("condition","race","bmi_raw"),
  condition_sex_bmi   = c("condition","newborn_sex","bmi_raw"),
  condition_race_sex   = c("condition","race","newborn_sex"),
  condition_race_sex_bmi = c("condition","race","newborn_sex","bmi_raw")
)

txi <- readRDS(here("data/processed_mrna/txi.rds"))
coldata <- readRDS(here("data/processed/metadata/metadata.rds"))

#designs <- list(condition_race_sex_bmi = c("condition","race","newborn_sex","bmi_raw"))
results_list_mRNA <- lapply(names(designs), function(name) {
  run_DE_analysis_mRNA(
    txi_obj = txi,
    colD = coldata,
    sample_prefix = paste0("mrna_cb_", name),
    design_variables = designs[[name]],
    alpha_v_ = 0.05,
    log2fc_thresh_v = 0.3,
    width = 2500,
    height = 2500
  )
})

names(results_list_mRNA) <- names(designs)

base::saveRDS(results_list_mRNA,here("data/processed_mrna/cb_results.rds"))

results_list_sig_mRNA <- lapply(names(results_list_mRNA),function(name) {
  get_sig(results_list_mRNA[[name]]$de$res,alpha = 0.05,lfc = 0.3)
})
names(results_list_sig_mRNA) <- names(designs)

cb_mrna_upset <- generate_upset_plt(results_list_sig_mRNA,"CB",type = "mRNA")
saveRDS(cb_mrna_upset,here("plots/cb_mrna_upset.rds"))
