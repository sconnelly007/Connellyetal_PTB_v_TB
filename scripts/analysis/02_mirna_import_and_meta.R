.libPaths("/nas/longleaf/home/sconnel/R/x86_64-pc-linux-gnu-library/4.4")
here::i_am("scripts/processing/02_process_mirna.Rmd")
library(tximeta)
library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(here)
library(DESeq2)
library(org.Hs.eg.db)
library(vroom)
set.seed(42)

samps_used <- vroom(here("data/raw/sample_sheet.csv"),delim = ",") %>% 
  mutate(pl_cb = if_else(str_detect(sample,"PL"),"placental","cord_blood"))

coldata <- readRDS(here("data/processed/metadata/metadata.rds")) %>% 
  filter(study_id != "044") %>% 
  filter(study_id != "017")

#cb
# remove 044, not in cord blood data
cb_coldata <- coldata %>% 
  filter(cb_match %in% samps_used$sample) %>% 
  mutate(cb_match = str_replace_all(cb_match,"-","."))

# Load counts 
mirna_cnts <- read_tsv(here("data/processed_mirna/mirna_quant/mirtop/mirna.tsv"))
colnames(mirna_cnts) <- str_replace(str_replace(colnames(mirna_cnts),"X",""),"_seqcluster","")
mirna_cnts <- mirna_cnts %>% column_to_rownames("miRNA")

#filter
mirna_cnts_CB <- mirna_cnts[,str_detect(colnames(mirna_cnts),"CB")]
# reorder to match columns
mirna_cnts_CB <- mirna_cnts_CB[,cb_coldata$cb_match]
saveRDS(list(counts = mirna_cnts_CB,coldata = cb_coldata),file = here("data/processed_mirna/cb_mirna.rds"))

#pl
pl_coldata <- coldata %>% 
  filter(pl_match %in% samps_used$sample) %>% 
  mutate(pl_match = str_replace_all(pl_match,"-","."))

#filter
mirna_cnts_PL <- mirna_cnts[,str_detect(colnames(mirna_cnts),"PL")]
# reorder to match columns
mirna_cnts_PL <- mirna_cnts_PL[,pl_coldata$pl_match]
saveRDS(list(counts = mirna_cnts_PL,coldata = pl_coldata),
        file = here("data/processed_mirna/pl_mirna.rds"))

