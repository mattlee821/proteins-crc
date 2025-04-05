rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(functions)
library(dplyr)
library(tidyr)
library(genetics.binaRies)
library(plinkbinr)
library(ieugwasr)
library(coloc)
library(parallel)
library(foreach)
library(doParallel)
source("/data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/functions_coloc-sensitivity.R")

# data ====
VAR_filelist <- "filelist-04"
LIST_FILES <- fread(paste0("/data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/filelist/", VAR_filelist), header = F, sep = "")
LIST_FILES <- LIST_FILES$V1

# VARS ====
priors <- list(
  list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-6), # 10k SNPs
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), # 1k SNPs
  list(p1 = 2e-5, p2 = 2e-5, p12 = 2e-6), # 5k SNPs
  list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-7) # conservative
)
label_priors <- c(
  "p1=1e-5;p2=1e-5;p12=1e-6", # 10k SNPs
  "p1=1e-4;p2=1e-4;p12=1e-5", # 1k SNPs
  "p1=2e-5;p2=2e-5;p12=2e-6", # 5k SNPs
  "p1=1e-5;p2=1e-5;p12=1e-7" # conservative
)

# df ====
data_coloc <- data.table::fread(LIST_FILES, header = T, sep = "\t", fill = TRUE)
data_coloc <- data_coloc %>%
  filter(!is.na(eaf.outcome) & !is.na(eaf.exposure))
# Skip if no significant p-values
if (!any(data_coloc$pval.exposure < 1E-5, na.rm = TRUE)) {
  return(NULL)
}
list_data <- split(x = data_coloc, f = data_coloc$id.outcome)

# data_coloc ====
tryCatch({
  if (!is.null(data_coloc)) {
    # reference ====
    VAR_population <- if (any(grepl("EUR", unique(data_coloc$id.exposure), ignore.case = TRUE))) "EUR" else "ALL"
    reference <- data.table::fread(paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/stats.stats"))
    
    # data ====
    data_coloc <- bind_rows(list_data) %>%
      rename(CHR = chr.exposure,
             POS = pos.exposure,
             EA = effect_allele.exposure,
             OA = other_allele.exposure,
             EAF = eaf.exposure,
             BETA = beta.exposure,
             SE = se.exposure,
             P = pval.exposure,
             N = samplesize.exposure,
             phenotype = id.exposure) %>%
      select(SNP, CHR, POS, EA, OA, EAF, BETA, SE, P, N, phenotype) %>%
      unique() %>%
      dplyr::inner_join(reference, by = c("SNP" = "Predictor")) %>%
      dplyr::mutate(BETA = ifelse(EA != A1, -BETA, BETA)) %>%  # Flip beta if EA != A1 to align GWAS with reference
      dplyr::distinct(SNP, .keep_all = TRUE)
    if (nrow(data_coloc) == 0) {
      return(NULL)
    }
    # ld ====
    ld <- ieugwasr::ld_matrix_local(
      data_coloc$SNP,
      with_alleles = FALSE, 
      bfile = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", 
                     VAR_population, "/", 
                     VAR_population),
      plink_bin = get_plink_exe())  
    if (is.null(ld) || nrow(ld) == 0) {
      return(NULL)
    }
    # data_coloc list ====
    data_coloc <- list(
      chr = data_coloc$CHR,
      position = data_coloc$POS,
      snp = rownames(ld),
      EA = data_coloc$EA,
      OA = data_coloc$OA,
      MAF = data_coloc$EAF,
      beta = data_coloc$BETA,
      SE = data_coloc$SE,
      pval = data_coloc$P,
      varbeta = data_coloc$SE^2,
      N = min(data_coloc$N),
      phenotype = unique(data_coloc$phenotype),
      Z = data_coloc$BETA / data_coloc$SE,
      LD = ld,
      type = "quant"
    )
    # end ====
  } else {
    message("list_data[[", i, "]] is NULL. Skipping this iteration.")
  }
}, error = function(e) {
  message("Error with alignment check and allele flipping for index ", i, ": ", e$message)
  return(NULL)
})
# Skip if no significant p-values
if (!any(data_coloc$pval < 1E-5, na.rm = TRUE)) {
  return(NULL)
}

# run parallel ====
## VARS_parallel ====
cl <- makeCluster(length(list_data))
registerDoParallel(cl)
clusterExport(cl, c("priors", "label_priors", "data_coloc"))
table_master <- data.frame()

# parallel processing ====
foreach(i = seq_along(list_data), .combine = rbind, .packages =
          c("data.table", "dplyr", "coloc", "functions", "genetics.binaRies", "plinkbinr")) %dopar% {
            # data_coloc2 ====
            tryCatch({
              if (!is.null(list_data[[i]])) {
                # reference ====
                VAR_population <- if (any(grepl("EUR", unique(list_data[[i]]$id.outcome), ignore.case = TRUE))) "EUR" else "ALL"
                reference <- data.table::fread(paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/stats.stats"))
                
                # data ====
                data_coloc2 <- list_data[[i]] %>%
                  filter(SNP %in% data_coloc$snp) %>% 
                  rename(CHR = chr.outcome,
                         POS = pos.exposure,
                         EA = effect_allele.outcome,
                         OA = other_allele.outcome,
                         EAF = eaf.outcome,
                         BETA = beta.outcome,
                         SE = se.outcome,
                         P = pval.outcome,
                         N = samplesize.outcome,
                         phenotype = id.outcome) %>%
                  select(SNP, CHR, POS, EA, OA, EAF, BETA, SE, P, N, phenotype) %>%
                  unique() %>%
                  dplyr::inner_join(reference, by = c("SNP" = "Predictor")) %>%
                  dplyr::mutate(BETA = ifelse(EA != A1, -BETA, BETA)) %>%  # Flip beta if EA != A1 to align GWAS with reference
                  dplyr::distinct(SNP, .keep_all = TRUE)
                if (nrow(data_coloc2) == 0) {
                  return(NULL)
                }
                # ld ====
                ld <- ieugwasr::ld_matrix_local(
                  data_coloc2$SNP,
                  with_alleles = FALSE, 
                  bfile = paste0("/data/GWAS_data/work/references/1000genomes/phase3/", 
                                 VAR_population, "/", 
                                 VAR_population),
                  plink_bin = get_plink_exe())  
                if (is.null(ld) || nrow(ld) == 0) {
                  return(NULL)
                }
                # data_coloc list ====
                data_coloc2 <- list(
                  chr = data_coloc2$CHR,
                  position = data_coloc2$POS,
                  snp = rownames(ld),
                  EA = data_coloc2$EA,
                  OA = data_coloc2$OA,
                  MAF = data_coloc2$EAF,
                  beta = data_coloc2$BETA,
                  SE = data_coloc2$SE,
                  pval = data_coloc2$P,
                  varbeta = data_coloc2$SE^2,
                  N = min(data_coloc2$N),
                  phenotype = unique(data_coloc2$phenotype),
                  Z = data_coloc2$BETA / data_coloc2$SE,
                  LD = ld,
                  type = "quant"
                )
              } else {
                message("list_data[[", i, "]] is NULL. Skipping this iteration.")
              }
            }, error = function(e) {
              message("Error with alignment check and allele flipping for index ", i, ": ", e$message)
              return(NULL)
            })
            
            # data_coloc1: match data_coloc to data_coloc2 ====
            snp_intersect <- intersect(data_coloc$snp, data_coloc2$snp)
            data_coloc1 <- lapply(names(data_coloc), function(name) {
              if (name == "snp") {
                # Filter SNPs to include only common SNPs
                data_coloc$snp[data_coloc$snp %in% snp_intersect]
              } else if (name == "LD") {
                # Subset LD matrix by SNPs (rows and columns)
                snp_indices <- which(data_coloc$snp %in% snp_intersect)
                data_coloc$LD[snp_indices, snp_indices, drop = FALSE]
              } else if (length(data_coloc[[name]]) == length(data_coloc$snp)) {
                # Subset other elements associated with each SNP
                data_coloc[[name]][data_coloc$snp %in% snp_intersect]
              } else {
                # Keep other elements (e.g., N, phenotype, type) unchanged
                data_coloc[[name]]
              }
            })
            data_coloc1 <- setNames(data_coloc1, names(data_coloc))
            data_coloc2$position <- data_coloc1$position # match pos for merging
            
            # finemap check ====
            SNP_causal_exposure <- coloc::finemap.abf(dataset = data_coloc1) %>%
              tidyr::drop_na() %>%
              filter(SNP.PP == max(SNP.PP)) %>%
              select(snp, SNP.PP)
            
            SNP_causal_outcome <- coloc::finemap.abf(dataset = data_coloc2) %>%
              tidyr::drop_na() %>%
              filter(SNP.PP == max(SNP.PP)) %>%
              select(snp, SNP.PP)
            
            # coloc.abf ====
            coloc_results <- lapply(priors, function(params) {
              coloc::coloc.abf(
                dataset1 = data_coloc1,
                dataset2 = data_coloc2,
                p1 = params$p1,
                p2 = params$p2,
                p12 = params$p12
              )
            })
            
            # results table ====
            table_coloc <- data.frame()
            # Loop through the coloc_results list
            for (k in seq_along(coloc_results)) {
              results_coloc <- coloc_results[[k]]
              results <- data.frame(
                id.exposure = data_coloc1$phenotype,
                id.outcome = data_coloc2$phenotype,
                finemap_snp_exposure = SNP_causal_exposure$snp,
                finemap_snp_exposure_PP = SNP_causal_exposure$SNP.PP,
                finemap_snp_outcome = SNP_causal_outcome$snp,
                finemap_snp_outcome_PP = SNP_causal_outcome$SNP.PP,
                nsnps = results_coloc$summary[["nsnps"]],
                h0 = results_coloc$summary[["PP.H0.abf"]],
                h1 = results_coloc$summary[["PP.H1.abf"]],
                h2 = results_coloc$summary[["PP.H2.abf"]],
                h3 = results_coloc$summary[["PP.H3.abf"]],
                h4 = results_coloc$summary[["PP.H4.abf"]],
                prior_p1 = results_coloc$priors[["p1"]],
                prior_p2 = results_coloc$priors[["p2"]],
                prior_p12 = results_coloc$priors[["p12"]],
                priors_label = label_priors[k]
              )
              # Bind the current results to the accumulated table
              table_coloc <- bind_rows(table_coloc, results)
            }
            # write ====
            directory <- paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/002_coloc/results/", data_coloc$phenotype, "/")
            if (!dir.exists(directory)) {
              dir.create(directory, 
                         recursive = TRUE)
            }
            write.table(
              table_coloc, 
              paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/002_coloc/results/", data_coloc$phenotype, "/", data_coloc2$phenotype, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
            
            # sensitivity  ====
            nsnps <- results_coloc$summary[["nsnps"]] # make plot based on number of SNPs
            id_nsnp <- apply(abs(outer(nsnps, c(10000, 1000, 5000), "-")), 1, which.min)
            id_nsnp <- c(id_nsnp, setdiff(1:length(priors), id_nsnp)) # Prioritize initial value
            plot_created <- FALSE  # Initialize a flag to track if the plot is created
            for (current_id in id_nsnp) {
              tryCatch({
                plot <- my_coloc_sensitivity(
                  obj = coloc_results[[current_id]], 
                  rule = "H4 > 0.8", 
                  npoints = 100, 
                  row = 1, 
                  suppress_messages = TRUE,
                  trait1_title = "trait 1", 
                  trait2_title = "trait 2",
                  dataset1 = NULL, 
                  dataset2 = NULL,
                  data_check_trait1 = data_coloc1, 
                  data_check_trait2 = data_coloc2
                )
                # If the plot is created successfully, update the flag and exit the inner loop
                plot_created <- TRUE
                id_nsnp <- current_id  # Update id_nsnp to the successful value
                break
              }, error = function(e) {
                # If an error occurs, do nothing and try the next id_nsnp
              })
            }
            # If no plot was created after trying all id_nsnp values, skip to the next iteration
            if (!plot_created) {
              message("Skipping to next item in the top-level loop.")
              return(NULL)
            }
            title <- cowplot::ggdraw() + 
              cowplot::draw_label(
                paste0(" trait 1 = ", data_coloc1$phenotype, "\n trait 2 = ", data_coloc2$phenotype, "\n priors = ", label_priors[id_nsnp]),
                fontface = 'bold',
                x = 0,
                hjust = 0
              )
            plot <- cowplot::plot_grid(
              title, plot,
              ncol = 1,
              rel_heights = c(0.1, 1)
            )
            
            ## save
            directory <- paste0("/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/002_coloc/figures/", data_coloc$phenotype, "/")
            if (!dir.exists(directory)) {
              dir.create(directory, 
                         recursive = TRUE)
            }
            tiff(filename = paste0(directory, data_coloc2$phenotype, ".tiff"), 
                 width = 1100, height = 500, units = "px")
            print(plot)
            dev.off()
            return(NULL)
            
          }

# Stop cluster
stopCluster(cl)
