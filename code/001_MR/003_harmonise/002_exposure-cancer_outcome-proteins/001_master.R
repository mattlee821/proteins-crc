rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(functions)
source("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/functions_harmonise.R")

VAR_filelist <- "filelist-01"
VAR_filelist <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/002_exposure-cancer_outcome-proteins/filelist/", VAR_filelist), header = F, sep = " ")
VAR_filelist <- VAR_filelist[[1]]
VAR_filelist <- basename(VAR_filelist)
VAR_filelist <- gsub(pattern = ".txt", replacement = "", x = VAR_filelist)
cat("# START \n")

# data ====
data_exposure <- fread("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/001_instruments/002_instruments-cancer/instruments.txt", header = T)
cat("## processing: ",VAR_filelist, "\n")
data_outcome <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/", VAR_filelist, ".txt"), header = T)

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
            file = paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/003_harmonise/002_exposure-cancer_outcome-proteins/", VAR_filelist, ".txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat("# FINISHED \n")
