rm(list = ls())
set.seed(821)

# environment ====
library(functions)
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
dyn.load("/home/leem/local/lib/libiconv.so.2")
Sys.setenv(LD_LIBRARY_PATH = "/home/leem/local/lib")
library(gwasvcf)
library(genetics.binaRies)
library(dplyr)
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
library(parallel)
library(foreach)
library(doParallel)
source("/data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/functions_harmonise.R")

# data: cancer ====
LIST_FILES <- c(
  list.files(path = "/data/GWAS_data/work/huyghe_2018_PMID30510241/GWAS/", pattern = ".gz", all.files = T, full.names = T),
  list.files(path = "/data/GWAS_data/work/laskar_2024_PMID38408508/GWAS/", pattern = ".gz", all.files = T, full.names = T),
  list.files(path = "/data/GWAS_data/work/fernandez-rozadilla_2022_PMID36539618/GWAS/", pattern = ".gz", all.files = T, full.names = T)
)
list_cancer <- lapply(X = LIST_FILES, FUN = data.table::fread, header = T, sep = "\t")

## VARS cancer ====
LABEL_STUDY <- sub(".*/([^/]+)/GWAS/.*\\.gz$", "\\1", LIST_FILES)
LABEL_COHORT <- ifelse(grepl("GECCO-UKB", LIST_FILES, ignore.case = TRUE), "GECCO-UKB",
                       ifelse(grepl("GECCO", LIST_FILES, ignore.case = TRUE), "GECCO", 
                              ifelse(grepl("UKB", LIST_FILES, ignore.case = TRUE), "UKB", 
                                     "GECCO")))
LABEL_SEX <- ifelse(grepl("combined", LIST_FILES, ignore.case = TRUE), "combined",
                    ifelse(grepl("female", LIST_FILES, ignore.case = TRUE), "female",
                           ifelse(grepl("male", LIST_FILES, ignore.case = TRUE), "male", NA)))
LABEL_PHENOTYPE <- ifelse(grepl("crc", LIST_FILES, ignore.case = TRUE), "overall",
                          ifelse(grepl("colon", LIST_FILES, ignore.case = TRUE), "colon",
                                 ifelse(grepl("rectal", LIST_FILES, ignore.case = TRUE), "rectal",
                                        ifelse(grepl("distal", LIST_FILES, ignore.case = TRUE), "distal",
                                               ifelse(grepl("proximal", LIST_FILES, ignore.case = TRUE), "proximal",
                                                      ifelse(grepl("early-onset", LIST_FILES, ignore.case = TRUE), "early-onset", NA))))))
LABEL_POPULATION <- ifelse(grepl("european", LIST_FILES, ignore.case = TRUE), "EUR",
                           ifelse(grepl("allancestries", LIST_FILES, ignore.case = TRUE), "ALL", 
                                  "ALL"))
VAR_cancer <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", LABEL_PHENOTYPE, ";", LABEL_COHORT)

# data: proteins ====
LIST_FILES <- c(
  # list.files(path = "/data/GWAS_data/work/ferkingstad_2021_PMID34857953/window/", pattern = ".txt", all.files = T, full.names = T),
  # list.files(path = "/data/GWAS_data/work/pietzner_2021_PMID34648354/window/", pattern = ".txt", all.files = T, full.names = T),
  # list.files(path = "/data/GWAS_data/work/sun_2023_PMID37794186/window/european/", pattern = ".txt", all.files = T, full.names = T),
  # list.files(path = "/data/GWAS_data/work/sun_2023_PMID37794186/window/combined/", pattern = ".txt", all.files = T, full.names = T),
  list.files(path = "/data/GWAS_data/work/zhang_2022_PMID35501419/GWAS/european/", pattern = ".txt", all.files = T, full.names = T)
)

## VARS proteins ====
LABEL_STUDY <- sub("^/data/GWAS_data/work/([^/]+)/.*$", "\\1", LIST_FILES)
LABEL_SEX <- "combined"
LABEL_COHORT <- ifelse(grepl("sun", LABEL_STUDY, ignore.case = TRUE), "UKB",
                       ifelse(grepl("ferkingstad", LABEL_STUDY, ignore.case = TRUE), "deCODE", 
                              ifelse(grepl("pietzner", LABEL_STUDY, ignore.case = TRUE), "FENLAND", 
                                     ifelse(grepl("zhang", LABEL_STUDY, ignore.case = TRUE), "ARIC", NA
                                     ))))
LABEL_POPULATION <- ifelse(grepl("european", LIST_FILES, ignore.case = TRUE), "EUR",
                           ifelse(grepl("combined", LIST_FILES, ignore.case = TRUE), "ALL", 
                                  "EUR"))
LABEL_PHENOTYPE <- basename(LIST_FILES)
LABEL_PHENOTYPE <- sub("\\..*", "", LABEL_PHENOTYPE)
VAR_protein <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", LABEL_PHENOTYPE, ";", LABEL_COHORT)

# run parallel ====
## VARS_parallel ====
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## parallel processing ====
foreach(i = seq_along(LIST_FILES), .packages = c("data.table", "dplyr", "TwoSampleMR", "functions")) %dopar% {
  # VARS ====
  listfile_name <- LIST_FILES[i]
  var_protein <- VAR_protein[i]
  
  # df_protein ====
  df_protein <- data.table::fread(listfile_name, header = T, sep = "\t", fill = TRUE)
  if (nrow(df_protein) == 0 || "X" %in% df_protein$CHR) {
    return(NULL)
    }
  
  # format df ====
  df_protein <- df_protein %>%
    as.data.frame() %>%
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
    )
  
  # VARS ====
  ID <- unique(df_protein$SNP)
  list_data <- list()
  
  # process df ====
  list_data <- lapply(list_cancer, function(df) {
    df %>%
      dplyr::filter(SNP %in% ID) %>%
      data.frame() %>%
      TwoSampleMR::format_data(
        dat = .,
        type = "outcome",
        snps = NULL,
        header = TRUE,
        phenotype_col = "phenotype",
        id_col = "phenotype",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        pval_col = "P",
        eaf_col = "EAF",
        effect_allele_col = "EA",
        other_allele_col = "OA",
        chr_col = "CHR",
        pos_col = "POS",
        samplesize_col = "N",
        min_pval = 1e-200,
        log_pval = FALSE
      )
  })

  # EAF ====
  list_data <- lapply(seq_along(list_data), function(j) {
    df_cancer <- list_data[[j]]
    
    # Check for missing EAF and replace missing
    if (anyNA(df_cancer$eaf.outcome)) {
      VAR_population <- if (any(grepl("EUR", VAR_cancer[j], ignore.case = TRUE))) "EUR" else "ALL"
      df_cancer <- functions::missing_EAF(
        df = df_cancer,
        reference = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/stats.stats"),
        column_EAF = "eaf.outcome"
      )
    }
    
    return(df_cancer)
  })    

  # harmonise ====
  list_data <- lapply(seq_along(list_data), function(j) {
    df_cancer <- list_data[[j]]
    
    # Harmonize data
    data_harmonise <- my_harmonise_data(
      exposure_dat = df_protein,
      outcome_dat = df_cancer,
      action = 2
    )
    
    # Add additional columns
    data_harmonise <- data_harmonise %>%
      mutate(
        exposure = var_protein,
        id.exposure = var_protein,
        outcome = VAR_cancer[j],
        id.outcome = VAR_cancer[j]
      )
    
    return(data_harmonise)
  })
  
  # write ====
  data <- dplyr::bind_rows(list_data)
  write.table(data, paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/001_window/", var_protein, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  return(NULL)
}

# Stop cluster
stopCluster(cl)
