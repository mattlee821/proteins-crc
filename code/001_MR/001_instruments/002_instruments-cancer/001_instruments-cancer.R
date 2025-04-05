rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(dplyr)
library(data.table)

# data ====
LIST_FILES <- c(
  list.files("/data/GWAS_data/work/huyghe_2018_PMID30510241/GWAS/", pattern = "gz", full.names = T),
  list.files("/data/GWAS_data/work/laskar_2024_PMID38408508/GWAS/", pattern = "gz", full.names = T),
  list.files("/data/GWAS_data/work/fernandez-rozadilla_2022_PMID36539618/GWAS/", pattern = "gz", full.names = T))

# extract ====
for (FILE in LIST_FILES){
  message(FILE)
  # extract sex from the file-name
  LABEL_STUDY <- sub(".*/([^/]+)/GWAS/.*\\.gz$", "\\1", FILE)
  LABEL_COHORT <- ifelse(grepl("GECCO-UKB", FILE, ignore.case = TRUE), "GECCO-UKB",
                         ifelse(grepl("GECCO", FILE, ignore.case = TRUE), "GECCO", 
                                ifelse(grepl("UKB", FILE, ignore.case = TRUE), "UKB",
                                "GECCO")))
  LABEL_SEX <- ifelse(grepl("combined", FILE, ignore.case = TRUE), "combined",
                      ifelse(grepl("female", FILE, ignore.case = TRUE), "female",
                             ifelse(grepl("male", FILE, ignore.case = TRUE), "male", NA)))
  LABEL_PHENOTYPE <- ifelse(grepl("crc", FILE, ignore.case = TRUE), "overall",
                         ifelse(grepl("colon", FILE, ignore.case = TRUE), "colon",
                                ifelse(grepl("rectal", FILE, ignore.case = TRUE), "rectal",
                                       ifelse(grepl("distal", FILE, ignore.case = TRUE), "distal",
                                              ifelse(grepl("proximal", FILE, ignore.case = TRUE), "proximal",
                                       ifelse(grepl("early-onset", FILE, ignore.case = TRUE), "early-onset", NA))))))
  LABEL_POPULATION <- ifelse(grepl("european", FILE, ignore.case = TRUE), "EUR",
                         ifelse(grepl("allancestries", FILE, ignore.case = TRUE), "ALL", 
                                "ALL"))
  LABEL_TEST <- "p_ld-0.001"
  LABEL <- paste0(LABEL_STUDY, ";", LABEL_POPULATION, ";", LABEL_SEX, ";", LABEL_PHENOTYPE, ";", LABEL_COHORT, ";", LABEL_TEST)
  
  # data read
  data <- fread(FILE, header = T, sep = "\t")
  data <- subset(data, P < 5E-08)
  
  # EAF: if EAF isn't present then extract it and add it
  if (!("EAF" %in% colnames(data))) {
    EAF <- fread(input = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", LABEL_POPULATION, "/stats.stats"),
                 header = TRUE,
                 select = c("Predictor", "A1", "A2", "MAF"), # Select only necessary columns
                 data.table = FALSE)
    EAF <- EAF %>%
      filter(Predictor %in% data$SNP) %>%
      rename(EAF = MAF)
    data <- left_join(data, EAF, by = c("SNP"="Predictor"))
    data <- data %>%
      mutate(EAF = ifelse(EA == A1 & OA == A2, EAF,   # If "EA"="A1" and "OA"="A2", keep EAF as it is
                          ifelse(EA == A2 & OA == A1, 1 - EAF, EAF))) %>%  # If "EA"="A2" and "OA"="A1", replace EAF with 1 - EAF
      select(-A1,-A2)
  }
  
  # Check if any NA values in the "EAF" column
  if (anyNA(data$EAF)) {
    data_eaf <- data[is.na(data$EAF), ]
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
        mutate(EAF = ifelse(EA == A1 & OA == A2, MAF,   
                            ifelse(EA == A2 & OA == A1, 1 - MAF, MAF))) %>%
        select(-A1, -A2, -MAF)
      # Remove rows from data based on SNP column in data_eaf
      data <- data[!data$SNP %in% data_eaf$SNP, ]
      # Merge data_eaf with data
      data <- rbind(data, data_eaf)
    }
  }
  
  ## format 
  data <- data %>%
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
    )
  
  # clump
  colnames(data)[colnames(data) == "SNP"] <- "rsid"
  colnames(data)[colnames(data) == "pval.exposure"] <- "pval"
  data <- ld_clump(dat = data,
                           clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                           pop = LABEL_POPULATION,
                           bfile = paste0("/data/GWAS_data/work/references/1000genomes/phase3/",LABEL_POPULATION,"/",LABEL_POPULATION),
                           plink_bin = genetics.binaRies::get_plink_binary())
  colnames(data)[colnames(data) == "rsid"] <- "SNP"
  colnames(data)[colnames(data) == "pval"] <- "pval.exposure"
  
  # if CHR/POS isn't present then extract it and add it
  if (!("chr.exposure" %in% colnames(data))) {
    info <- fread(input = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", LABEL_POPULATION, "/", LABEL_POPULATION, ".bim"),
                 header = FALSE, col.names = c("chr.exposure", "SNP", "V1", "pos.exposure", "V2", "V3"))
    info <- info %>%
      filter(SNP %in% data$SNP) %>%
      select("SNP", "chr.exposure", "pos.exposure")
    data <- left_join(data, info, by = c("SNP"="SNP"))
  }
  data$chr.exposure <- as.integer(data$chr.exposure)
  
  ## save
  data$id.exposure <- LABEL
  data$exposure <- LABEL
  write.table(data, paste0("analysis/001_MR/001_instruments/002_instruments-cancer/", LABEL, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# combine ====
list_data <- list.files(path = "analysis/001_MR/001_instruments/002_instruments-cancer/", pattern = "txt", all.files = T, full.names = T)
list_data <- lapply(list_data, fread, header = T, sep = "\t")
data <- bind_rows(list_data)
write.table(data, "analysis/001_MR/001_instruments/002_instruments-cancer/instruments.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
