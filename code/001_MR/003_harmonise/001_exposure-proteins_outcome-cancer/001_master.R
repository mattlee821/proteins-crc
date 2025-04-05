rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(functions)
source("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/functions_harmonise.R")

VAR_filelist <- "filelist-01"
VAR_filelist <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/001_exposure-proteins_outcome-cancer/filelist/", VAR_filelist), header = F, sep = " ")
VAR_filelist <- VAR_filelist[[1]]
cat("# START \n")

# data ====
data_exposure <- fread("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/001_instruments/001_instruments-proteins/instruments.txt", header = T)
cat("## processing: ",VAR_filelist, "\n")
data_exposure <- data_exposure[data_exposure$id.exposure == VAR_filelist, ]
cat("## N SNPs: ", nrow(data_exposure), "\n")
data_outcome <- list.files(path = "/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/001_exposure-proteins_outcome-cancer/", pattern = ".txt", all.files = T, full.names = T)
data_outcome <- lapply(data_outcome, fread, header = TRUE, sep = "\t")
data_outcome <- bind_rows(data_outcome)

# harmonise ====
cat("## harmonising \n")
data_harmonise <- my_harmonise_data( # this uses lef_join instead of merge
  exposure_dat = data_exposure, 
  outcome_dat = data_outcome, 
  action = 2)
cat("## harmonised \n")
data_harmonise <- data_harmonise %>%
  select(-pval_origin.outcome, -data_source.outcome, -pval_origin.exposure)

# write ====
cat("## write \n")
write.table(x = data_harmonise, 
            file = paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/003_harmonise/001_exposure-proteins_outcome-cancer/", VAR_filelist, ".txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat("# FINISHED \n")
