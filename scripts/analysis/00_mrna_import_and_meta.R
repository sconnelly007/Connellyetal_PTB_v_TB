.libPaths("/nas/longleaf/home/sconnel/R/x86_64-pc-linux-gnu-library/4.4")
here::i_am("scripts/processing/01_process_mrna.Rmd")
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

samps_used_mrna <- vroom(here("data/raw/sample_sheet_mRNA.csv"),delim = ",")

# import more metadata
coldata <- vroom(here("data","metadata","VoraPretermBirthAndP_DATA_LABELS_2023-06-29_1205.csv")) %>% 
  janitor::clean_names()
coldata <- coldata %>% 
  dplyr::select(study_id,study_arm,race,ga_at_delivery_257,ga_at_delivery_258,best_edc,delivery_date_time,preterm_delivery_delivery_ga_37_weeks,indication_for_preterm_delivery,newborn_sex,bmi) %>% 
  filter(!(study_id == "077")) %>% 
  mutate(study_id = case_when(
    study_id == "33" ~ "033",
    study_id == "54" ~ "054",
    study_id == "77" ~ "077",
    study_id == "79" ~ "079",
    study_id == "87" ~ "087",
    study_id == "99" ~ "099",
    study_id == "BIG 012" ~ "BG12",
    .default = as.character(study_id)
  )) %>% 
  rename("ga_at_delivery_weeks" = ga_at_delivery_257,
         "ga_at_delivery_days" = ga_at_delivery_258) %>% 
  mutate(is_ptb = if_else(ga_at_delivery_weeks < 37, TRUE,FALSE),
         race = str_replace(race," ","_")) %>% 
  mutate(bmi_raw = bmi)

# append 13-CB onto them
coldata <- coldata %>% 
  mutate(cb_match = paste0("13-CB-",study_id),
         pl_match = paste0("13-PL-",study_id))

# you need to filter out the sample missing BMI
coldata <- coldata %>% 
  filter(cb_match %in% samps_used_mrna$sample)
dir <- here("data", "processed_mrna")
coldata$files <- file.path(dir, "star_salmon", coldata$cb_match, "quant.sf")
coldata$cb_match <- coldata$cb_match %>% 
  factor()
coldata$pl_match <- coldata$pl_match %>% 
  factor()
coldata$condition <- coldata$is_ptb %>%
  factor(levels = c(FALSE,TRUE)) %>%
  fct_recode( "Term" = "FALSE","Preterm" = "TRUE")
coldata$names <- coldata$cb_match
coldata$newborn_sex <- as.factor(coldata$newborn_sex)
coldata <- coldata %>% 
  mutate(ga_at_delivery_combined = ga_at_delivery_weeks + ga_at_delivery_days / 7) %>% 
  mutate(ga_bin = case_when(
    ga_at_delivery_combined < 28 ~ "<28w (extreme preterm)",
    ga_at_delivery_combined < 32 ~ "28–31w (very preterm)",
    ga_at_delivery_combined < 34 ~ "32–33w (moderate preterm)",
    ga_at_delivery_combined < 37 ~ "34–36w (late preterm)",
    ga_at_delivery_combined < 42 ~ "37–41w (term)",
    TRUE ~ NA_character_
  )) %>% 
  mutate(bmi = case_when(
    bmi_raw > 19 & bmi_raw < 30 ~ "19-30",
    bmi_raw > 30 & bmi_raw < 40 ~ "30-40",
    bmi_raw > 40 & bmi_raw < 50 ~ "40-50",
    bmi_raw > 50 ~ ">50",
    TRUE ~ NA_character_
  ))
coldata$bmi <- coldata$bmi %>% 
  factor(levels = c("19-30","30-40","40-50",">50"))
coldata$ga_bin <- coldata$ga_bin %>% 
  factor(levels = c("<28w (extreme preterm)","28–31w (very preterm)","32–33w (moderate preterm)",
                    "34–36w (late preterm)","37–41w (term)"))
coldata$race <- as.factor(coldata$race)
coldata$bmi_raw <- scale(coldata$bmi_raw, center=TRUE, scale=TRUE)
coldata <- coldata %>% 
  mutate(indication_for_preterm_delivery = ifelse(grepl(pattern = "PPROM",x = indication_for_preterm_delivery),
                                                  "PPROM",
                                                  ifelse(grepl(pattern = "PTL",x = indication_for_preterm_delivery),
                                                         "PTL","Term")))

coldata <- coldata %>% 
  filter(cb_match != "13-CB-044") %>% 
  # case, unfortunately
  # 13-CB-017
  filter(cb_match != "13-CB-017")

rownames(coldata) <- coldata$names

saveRDS(coldata,
        here("data/processed/metadata/metadata.rds"))

tx2gene <- read_tsv(
  here("data", "processed_mrna", "star_salmon", "salmon.merged.tx2gene.tsv")
) %>% 
  select(transcript_id,gene_id)

txi <- tximport::tximport(files = coldata$files,
                         type = "salmon",
                         tx2gene = tx2gene)

saveRDS(txi,
        here("data/processed_mrna/txi.rds"))
