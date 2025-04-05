rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(data.table)
# Define a function to load and process the data
process_data <- function(VAR, pattern, path_prefix) {
  list_files <- list.files(
    paste0("/data/GWAS_data/work/", VAR, path_prefix), 
    pattern = pattern, all.files = TRUE, full.names = TRUE
  )
  
  list_data <- lapply(list_files, fread, header = TRUE, sep = "\t")
  
  data <- bind_rows(list_data) %>%
    mutate(
      EA = toupper(EA),   # Convert EA to uppercase
      OA = toupper(OA)    # Convert OA to uppercase
    ) %>%
    filter(
      BETA != 0,
      EA %in% c("A", "C", "T", "G"),
      OA %in% c("A", "C", "T", "G") 
    ) %>%
    mutate(
      ID = ifelse(
        test == "p_ld", 
        paste(VAR, phenotype, paste(test, LD, sep = "-"), sep = ";"), 
        paste(VAR, phenotype, test, sep = ";")
      )
    ) %>%
    mutate(
      Plog10 = ifelse(
        !"Plog10" %in% colnames(data), 
        -log10(P),  # Create Plog10 column if it doesn't exist
        ifelse(is.na(Plog10) | is.infinite(Plog10), -log10(P), Plog10)  # Handle NA/Inf values
      )
    )
  
  # Filter the SuSiE and Finimom credible sets
  temp <- data %>%
    filter(test %in% c("SuSiE", "Finimom")) %>%
    group_by(ID, cs) %>%
    filter(
      (test == "Finimom" & !all(is.na(P)) & P == min(P, na.rm = TRUE)) |
        (test == "SuSiE" & !all(is.na(PIP)) & PIP == min(PIP, na.rm = TRUE) & 
           P == min(P[PIP == min(PIP, na.rm = TRUE)], na.rm = TRUE))
    ) %>%
    slice(1) %>%
    ungroup()
  
  data <- data %>%
    filter(!test %in% c("SuSiE", "Finimom")) %>%
    bind_rows(temp)
  cat("#", VAR, "\n")
  cat("# N phenotype:", length(unique(data$phenotype)), "\n")
  cat("# N ID:", length(unique(data$ID)), "\n")
  
  return(data)
}

# ferkingstad_2021_PMID34857953 ====
LABEL_STUDY <- "ferkingstad_2021_PMID34857953"
LABEL_SEX <- "combined"
LABEL_POPULATION <- "EUR"
LABEL_COHORT <- "deCODE"
data_ferkingstad <- process_data(VAR = LABEL_STUDY, pattern = ".txt", path_prefix = "/finemapping/data/")
data_ferkingstad$phenotype <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", data_ferkingstad$phenotype, ";", LABEL_COHORT, ";", sapply(strsplit(data_ferkingstad$ID, ";"), function(x) tail(x, 1)))

# pietzner_2021_PMID34648354 ====
LABEL_STUDY <- "pietzner_2021_PMID34648354"
LABEL_SEX <- "combined"
LABEL_POPULATION <- "EUR"
LABEL_COHORT <- "FENLAND"
data_pietzner <- process_data(VAR = "pietzner_2021_PMID34648354", pattern = ".txt", path_prefix = "/finemapping/data/")
data_pietzner$phenotype <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", data_pietzner$phenotype, ";", LABEL_COHORT, ";", sapply(strsplit(data_pietzner$ID, ";"), function(x) tail(x, 1)))

# sun_2023_PMID37794186: EUR ====
LABEL_STUDY <- "sun_2023_PMID37794186"
LABEL_SEX <- "combined"
LABEL_POPULATION <- "EUR"
LABEL_COHORT <- "UKB"
data_sun_EUR <- process_data(VAR = "sun_2023_PMID37794186", pattern = ".txt", path_prefix = "/finemapping/european/data/")
data_sun_EUR$phenotype <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", data_sun_EUR$phenotype, ";", LABEL_COHORT, ";", sapply(strsplit(data_sun_EUR$ID, ";"), function(x) tail(x, 1)))

# sun_2023_PMID37794186: ALL ====
LABEL_STUDY <- "sun_2023_PMID37794186"
LABEL_SEX <- "combined"
LABEL_POPULATION <- "ALL"
LABEL_COHORT <- "UKB"
data_sun_ALL <- process_data(VAR = "sun_2023_PMID37794186", pattern = ".txt", path_prefix = "/finemapping/combined/data/")
data_sun_ALL$phenotype <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", data_sun_ALL$phenotype, ";", LABEL_COHORT, ";", sapply(strsplit(data_sun_ALL$ID, ";"), function(x) tail(x, 1)))

# zhang_2022_PMID35501419 ====
LABEL_STUDY <- "zhang_2022_PMID35501419"
LABEL_SEX <- "combined"
LABEL_POPULATION <- "EUR"
LABEL_COHORT <- "ARIC"
data_zhang <- process_data(VAR = "zhang_2022_PMID35501419", pattern = ".txt", path_prefix = "/finemapping/european/data/")
data_zhang$phenotype <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", data_zhang$phenotype, ";", LABEL_COHORT, ";", sapply(strsplit(data_zhang$ID, ";"), function(x) tail(x, 1)))

# combine ====
data <- bind_rows(data_ferkingstad, data_pietzner, data_sun_EUR, data_sun_ALL, data_zhang) %>%
  data.frame() %>%
  TwoSampleMR::format_data(
    dat = ., 
    type = "exposure", 
    phenotype_col = "phenotype", 
    snp_col = "SNP", 
    beta_col = "BETA", 
    se_col = "SE", 
    eaf_col = "EAF", 
    effect_allele_col = "EA", 
    other_allele_col = "OA", 
    pval_col = "P", 
    samplesize_col = "N", 
    chr_col = "CHR", 
    pos_col = "POS38"
  ) %>%
  mutate(id.exposure = exposure)

write.table(data, "analysis/001_MR/001_instruments/001_instruments-proteins/instruments.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# VEP ====
temp <- data %>%
  dplyr::select(chr.exposure, pos.exposure, effect_allele.exposure, other_allele.exposure) %>%
  dplyr::mutate(driver = 1) %>%
  dplyr::mutate(chr.exposure = paste0("chr", chr.exposure)) %>%
  unique()
write.table(temp, "analysis/001_MR/001_instruments/001_instruments-proteins/variants_with_driver_stat.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
