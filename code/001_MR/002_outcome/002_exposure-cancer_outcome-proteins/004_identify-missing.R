rm(list=ls())
set.seed(821)

# environment ====
dyn.load("/home/leem/local/lib/libiconv.so.2")
Sys.setenv(LD_LIBRARY_PATH = "/home/leem/local/lib")
library(functions)
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(gwasvcf)
library(genetics.binaRies)
library(dplyr)
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

# filelists ====
VAR_filelist <- list.files(
  path = "/data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/filelist/", 
  pattern = "^filelist-[0-9]+$", 
  all.files = TRUE, 
  full.names = TRUE
)

VAR_FILENAME <- c()
VAR_LABEL <- c()
formatted_count <- 0
# Iterate through each filelist
for (VAR_file in VAR_filelist) {
  formatted_count <- formatted_count + 1
  cat("Formatting filelist:", formatted_count, "of", length(VAR_filelist), "\n")
  
  # Read the list of files from the filelist
  list_files <- data.table::fread(
    VAR_file, 
    header = FALSE
  )[[1]]
  
  # Process each file in the list
  for (FILE in list_files) {
    # Extract label components
    LABEL_STUDY <- sub(".*/([^/]+)/GWAS/.*\\.gz$", "\\1", FILE)
    LABEL_SEX <- "combined"
    LABEL_COHORT <- ifelse(grepl("sun", LABEL_STUDY, ignore.case = TRUE), "UKB",
                           ifelse(grepl("ferkingstad", LABEL_STUDY, ignore.case = TRUE), "deCODE", 
                                  ifelse(grepl("pietzner", LABEL_STUDY, ignore.case = TRUE), "FENLAND", 
                                         ifelse(grepl("zhang", LABEL_STUDY, ignore.case = TRUE), "ARIC", NA
                                         ))))
    LABEL_POPULATION <- ifelse(grepl("european", FILE, ignore.case = TRUE), "EUR",
                               ifelse(grepl("combined", FILE, ignore.case = TRUE), "ALL", 
                                      "EUR"))
    
    # Read a subset of the file to extract phenotype
    data <- readLines(FILE, n = 2)
    header <- strsplit(data[1], "\t")[[1]]
    values <- strsplit(data[2], "\t")[[1]]
    # Create the dataframe
    data <- as.data.frame(matrix(values, nrow = 1, dimnames = list(NULL, header)), stringsAsFactors = FALSE)

    # Construct VAR_FILENAME and VAR_LABEL
    VAR_FILENAME <- c(VAR_FILENAME, FILE)
    VAR_LABEL <- c(VAR_LABEL, paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", unique(data$phenotype), ";", LABEL_COHORT))
  }
}

# created files ====
VAR_files <- list.files(
  path = "/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/", 
  pattern = "txt", 
  all.files = TRUE
) 
VAR_files <- gsub("\\.txt$", "", VAR_files)

# missing ====
VAR_setdiff <- setdiff(VAR_LABEL, VAR_files)

# make file path name ====
VAR_file <- c()
# Loop through each value in VAR_setdiff
for (i in seq_along(VAR_setdiff)) {
  # Split the current value by ";"
  parts <- strsplit(VAR_setdiff[i], ";")[[1]]
  
  # Logic to determine the file path
  if (grepl("ferkingstad", parts[1])) {
    file_path <- paste0("/data/GWAS_data/work/", parts[1], "/GWAS/", parts[4], ".txt.gz")
  } else if (grepl("sun", parts[1]) && parts[2] == "ALL") {
    file_path <- paste0("/data/GWAS_data/work/", parts[1], "/GWAS/combined/", parts[4], ".gz")
  } else if (grepl("sun", parts[1]) && parts[2] == "EUR") {
    file_path <- paste0("/data/GWAS_data/work/", parts[1], "/GWAS/european/", parts[4], ".gz")
  } else {
    file_path <- NA  # Default case if no condition is met
  }
  
  # Append the file path to VAR_file
  VAR_file <- c(VAR_file, file_path)
}

# process missing ====
cat("# STARTED \n")
# data ====
data_exposure <- data.table::fread("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/001_instruments/002_instruments-cancer/instruments.txt", header = T, sep = "\t")

# extarct data ====
for (FILE in VAR_file){
  LABEL_STUDY <- sub(".*/([^/]+)/GWAS/.*\\.gz$", "\\1", FILE)
  LABEL_SEX <- "combined"
  LABEL_COHORT <- ifelse(grepl("sun", LABEL_STUDY, ignore.case = TRUE), "UKB",
                         ifelse(grepl("ferkingstad", LABEL_STUDY, ignore.case = TRUE), "deCODE", 
                                ifelse(grepl("pietzner", LABEL_STUDY, ignore.case = TRUE), "FENLAND", 
                                       ifelse(grepl("zhang", LABEL_STUDY, ignore.case = TRUE), "ARIC", NA
                                       ))))
  LABEL_POPULATION <- ifelse(grepl("european", FILE, ignore.case = TRUE), "EUR",
                             ifelse(grepl("combined", FILE, ignore.case = TRUE), "ALL", 
                                    "EUR"))
  
  # data read
  data <- fread(FILE, header = T, sep = "\t") 
  data <- data[data$SNP %in% data_exposure$SNP]
  LABEL <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", unique(data$phenotype), ";", LABEL_COHORT)
  data$phenotype <- LABEL
  
  # data format
  data <- as.data.frame(data)
  data <- format_data(
    data,
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
  data <- data %>%
    filter(mr_keep.outcome == TRUE)
  
  # proxy search
  data <- functions::proxy_search(data_exposure = data_exposure, 
                                  data_outcome = data, 
                                  data_outcome_path = FILE, 
                                  outcome_sep = "\t",
                                  outcome_SNP = "SNP",
                                  outcome_BETA = "BETA",
                                  outcome_SE = "SE",
                                  outcome_P = "P",
                                  outcome_EA = "EA",
                                  outcome_OA = "OA",
                                  outcome_EAF = "EAF",
                                  outcome_N = "N",
                                  outcome_CHR = "CHR",
                                  outcome_POS = "POS38",
                                  outcome_ID = "phenotype",
                                  outcome_phenotype = "phenotype",
                                  data_reference = paste0("/data/GWAS_data/work/references/1000genomes/phase3/",LABEL_POPULATION,"/",LABEL_POPULATION,".bim"), 
                                  data_reference_path = paste0("/data/GWAS_data/work/references/1000genomes/phase3/",LABEL_POPULATION,"/",LABEL_POPULATION),
                                  tag_r2 = 0.8, 
                                  tag_kb = 5000, 
                                  tag_nsnp = 5000)
  data <- data[!duplicated(data$SNP), ]
  
  # EAF: if EAF isn't present then extract it and add it
  if (!any(c("EAF", "eaf.outcome") %in% colnames(data))) {
    EAF <- fread(input = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", LABEL_POPULATION, "/stats.stats"),
                 header = TRUE,
                 select = c("Predictor", "A1", "A2", "MAF"), # Select only necessary columns
                 data.table = FALSE)
    EAF <- EAF %>%
      filter(Predictor %in% data$SNP) %>%
      rename(EAF = MAF)
    data <- left_join(data, EAF, by = c("SNP"="Predictor"))
    data <- data %>%
      mutate(EAF = ifelse(effect_allele.outcome == A1 & other_allele.outcome == A2, EAF,   # If "EA"="A1" and "OA"="A2", keep EAF as it is
                          ifelse(effect_allele.outcome == A2 & other_allele.outcome == A1, 1 - EAF, EAF))) %>%  # If "EA"="A2" and "OA"="A1", replace EAF with 1 - EAF
      select(-A1,-A2) %>%
      rename(eaf.outcome = EAF)
  }
  # Check if any NA values in the "EAF" column
  if (anyNA(data$eaf.outcome)) {
    data_eaf <- data[is.na(data$eaf.outcome), ]
    # Apply your code block to update the data frame
    EAF <- fread(input = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", LABEL_POPULATION, "/stats.stats"),
                 header = TRUE,
                 select = c("Predictor", "A1", "A2", "MAF"), # Select only necessary columns
                 data.table = FALSE)
    # Filter EAF data
    EAF <- EAF %>%
      filter(Predictor %in% data_eaf$SNP)
    # Check if EAF data frame has any rows
    if (nrow(EAF) > 0) {
      # Calculate EAF
      data_eaf <- left_join(data_eaf, EAF, by = c("SNP" = "Predictor")) %>%
        mutate(eaf.outcome = ifelse(effect_allele.outcome == A1 & other_allele.outcome == A2, MAF,   
                                    ifelse(effect_allele.outcome == A2 & other_allele.outcome == A1, 1 - MAF, MAF))) %>%
        select(-A1, -A2, -MAF)
      # Remove rows from data based on SNP column in data_eaf
      data <- data[!data$SNP %in% data_eaf$SNP, ]
      # Merge data_eaf with data
      data <- rbind(data, data_eaf)
    }
  }
  
  # save
  data$id.outcome <- LABEL
  data$outcome <- LABEL
  write.table(data, 
              paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/", LABEL, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

cat("# FINISHED \n")



