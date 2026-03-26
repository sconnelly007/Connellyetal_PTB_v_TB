designs <- list(
  condition        = c("condition"),
  condition_race   = c("condition","race"),
  condition_sex    = c("condition","newborn_sex"),
  condition_bmi   = c("condition","bmi"),
  condition_race_bmi   = c("condition","race","bmi"),
  condition_sex_bmi   = c("condition","newborn_sex","bmi"),
  condition_race_sex   = c("condition","race","newborn_sex"),
  condition_race_sex_bmi = c("condition","race","newborn_sex","bmi")
)

# CB miRNA ----------------------------------------------------------------
cb_mirna <- readRDS(here("data/processed_mirna/cb_mirna.rds"))

results_list_CB <- lapply(names(designs), function(name) {
  run_DE_analysis_mirna(
    counts = cb_mirna$counts,
    coldata = cb_mirna$coldata,
    sample_prefix = paste0("mirna_mirtop_cb_", name),
    design_variables = designs[[name]],
    alpha_v_ = 0.05,
    log2fc_thresh_v = 0.3
  )
})
names(results_list_CB) <- names(designs)

base::saveRDS(results_list_CB,
        here("data/processed_mirna/cb_mirna_results.rds"))

results_list_sig_CB <- lapply(names(results_list_CB),function(name) {
  get_sig(results_list_CB[[name]]$de$res,alpha = 0.05,lfc = 0.3)
})
names(results_list_sig_CB) <- names(designs)

cb_upset <- generate_upset_plt(results_list_sig_CB,"CB",type = "miRNA")
saveRDS(cb_upset,here("plots/cb_mirna_upset.rds"))
# PL miRNA ----------------------------------------------------------------
pl_mirna <- readRDS(here("data/processed_mirna/pl_mirna.rds"))

results_list_pl <- lapply(names(designs), function(name) {
  run_DE_analysis_mirna(
    counts = pl_mirna$counts,
    coldata = pl_mirna$coldata,
    sample_prefix = paste0("mirna_mirtop_pl_", name),
    design_variables = designs[[name]],
    alpha_v_ = 0.05,
    log2fc_thresh_v = 0.3
  )
})

names(results_list_pl) <- names(designs)

base::saveRDS(results_list_pl,
        here("data/processed_mirna/pl_results.rds"))

results_list_sig_pl <- lapply(names(results_list_pl),function(name) {
  get_sig(results_list_pl[[name]]$de$res,alpha = 0.05,lfc = 0.3)
})
names(results_list_sig_pl) <- names(designs)

pl_upset <- generate_upset_plt(results_list_sig_pl,"PL",type = "miRNA")
saveRDS(pl_upset,here("plots/pl_mirna_upset.rds"))
# any shared? -------------------------------------------------------------
intersect(results_list_sig_CB$condition_race_sex_bmi, results_list_sig_pl$condition_race_sex_bmi)

intersect_df <- list(`CB` = results_list_sig_CB$condition_race_sex_bmi, 
                     `PL` = results_list_sig_pl$condition_race_sex_bmi)

generate_upset_plt(intersect_df,"shared","miRNA")
