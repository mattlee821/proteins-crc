rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(functions)

# source ====
cat("# START \n")
VAR_filelist <- "filelist-01"
VAR_filelist <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/004_MR/002_exposure-cancer_outcome-proteins/filelist/", VAR_filelist), header = F, sep = " ")
VAR_filelist <- VAR_filelist[[1]]
VAR_filelist <- basename(VAR_filelist)
VAR_filelist <- gsub(pattern = ".txt", replacement = "", x = VAR_filelist)

cat("## processing: ",VAR_filelist, "\n")

# data ====
cat("## data \n")
data_harmonise <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/003_harmonise/002_exposure-cancer_outcome-proteins/", VAR_filelist, ".txt"))

# data VARS ====
VAR_outcome_samplesize <- data_harmonise %>%
  select(id.outcome, samplesize.outcome) %>%
  unique() %>%
  summarise(min_samplesize = min(samplesize.outcome, na.rm = TRUE)) %>%
  pull(min_samplesize)

VAR_exposure_samplesize <- data_harmonise %>%
  select(id.exposure, samplesize.exposure) %>%
  unique()

# MR VARS ====
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods_heterogeneity <- setdiff(methods_heterogeneity, c("mr_ivw_radial"))
methods <- methods$obj
methods <- setdiff(methods, c("mr_raps", "mr_ivw_radial"))

# MR ====
cat("## MR \n")
data_mr <- mr(dat = data_harmonise, method_list = methods)
## format
data_mr <- data_mr %>%
  separate(
    exposure, 
    into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), 
    sep = ";", remove = TRUE) %>%
  separate(
    outcome, 
    into = c("outcome_study", "outcome_population", "outcome_sex", "outcome", "outcome_cohort"), 
    sep = ";", remove = TRUE) %>%
  left_join(., VAR_exposure_samplesize, by = "id.exposure") %>%
  rename(exposure_samplesize = samplesize.exposure) %>%
  mutate(outcome_samplesize = VAR_outcome_samplesize) %>%
  select(
    id.exposure, id.outcome,
    exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
    outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize,
    method, nsnp, b, se, pval
  ) %>%
  rename(
    BETA = b,
    SE = se,
    P = pval
  )
## save
write.table(
  data_mr, 
  paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/MR/", VAR_filelist, ".txt"), 
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# single snp ====
cat("## single-snp \n")
data_sensitivity <- mr_singlesnp(data_harmonise)
## format
data_sensitivity <- data_sensitivity %>%
  filter(!str_detect(SNP, regex("All", ignore_case = TRUE))) %>%
  separate(
    exposure, 
    into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), 
    sep = ";", remove = TRUE) %>%
  separate(
    outcome, 
    into = c("outcome_study", "outcome_population", "outcome_sex", "outcome", "outcome_cohort"), 
    sep = ";", remove = TRUE) %>%
  left_join(., VAR_exposure_samplesize, by = "id.exposure") %>%
  rename(exposure_samplesize = samplesize.exposure) %>%
  mutate(outcome_samplesize = VAR_outcome_samplesize) %>%
  select(
    id.exposure, id.outcome,
    exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
    outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize,
    SNP, b, se, p
  ) %>%
  rename(
    BETA = b,
    SE = se,
    P = p
  )
## save
write.table(
  data_sensitivity, 
  paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/singlesnp/", VAR_filelist, ".txt"), 
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# hetrogeneity ====
if (any(data_mr$nsnp > 3)) {
  cat("## hetrogeneity \n")
  data_sensitivity <- mr_heterogeneity(data_harmonise, method_list = methods_heterogeneity)
  ## format
  data_sensitivity <- data_sensitivity %>%
    separate(
      exposure, 
      into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), 
      sep = ";", remove = TRUE) %>%
    separate(
      outcome, 
      into = c("outcome_study", "outcome_population", "outcome_sex", "outcome", "outcome_cohort"), 
      sep = ";", remove = TRUE) %>%
    left_join(., VAR_exposure_samplesize, by = "id.exposure") %>%
    rename(exposure_samplesize = samplesize.exposure) %>%
    mutate(outcome_samplesize = VAR_outcome_samplesize) %>%
    select(
      id.exposure, id.outcome,
      exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
      outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize,
      method, Q, Q_df, Q_pval
    ) %>%
    rename(
      Q_P = Q_pval
    )
  ## save
  write.table(
    data_sensitivity, 
    paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/heterogeneity/", VAR_filelist, ".txt"), 
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# pleiotropy ====
if (any(data_mr$nsnp > 3)) {
  cat("## pleiotropy \n")
  data_sensitivity <- mr_pleiotropy_test(data_harmonise)
  ## format
  data_sensitivity <- data_sensitivity %>%
    separate(
      exposure, 
      into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), 
      sep = ";", remove = TRUE) %>%
    separate(
      outcome, 
      into = c("outcome_study", "outcome_population", "outcome_sex", "outcome", "outcome_cohort"), 
      sep = ";", remove = TRUE) %>%
    left_join(., VAR_exposure_samplesize, by = "id.exposure") %>%
    rename(exposure_samplesize = samplesize.exposure) %>%
    mutate(outcome_samplesize = VAR_outcome_samplesize) %>%
    select(
      id.exposure, id.outcome,
      exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
      outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize,
      egger_intercept, se, pval
    ) %>%
    rename(
      SE = se,
      P = pval
    )
  ## save
  write.table(
    data_sensitivity, 
    paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/pleiotropy/", VAR_filelist, ".txt"), 
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# leave one out ====
if (any(data_mr$nsnp > 3)) {
  cat("## leae-one-out \n")
  data_sensitivity <- mr_leaveoneout(data_harmonise)
  ## format
  data_sensitivity <- data_sensitivity %>%
    separate(
      exposure, 
      into = c("exposure_study", "exposure_population", "exposure_sex", "exposure", "exposure_cohort", "exposure_finemap"), 
      sep = ";", remove = TRUE) %>%
    separate(
      outcome, 
      into = c("outcome_study", "outcome_population", "outcome_sex", "outcome", "outcome_cohort"), 
      sep = ";", remove = TRUE) %>%
    left_join(., VAR_exposure_samplesize, by = "id.exposure") %>%
    rename(exposure_samplesize = samplesize.exposure) %>%
    mutate(outcome_samplesize = VAR_outcome_samplesize) %>%
    select(
      id.exposure, id.outcome,
      exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
      outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize,
      SNP, b, se, p
    ) %>%
    rename(
      BETA = b,
      SE = se,
      P = p
    )
  ## save
  write.table(
    data_sensitivity, 
    paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/leaveoneout/", VAR_filelist, ".txt"), 
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# END ====
cat("# FINISHED \n")
