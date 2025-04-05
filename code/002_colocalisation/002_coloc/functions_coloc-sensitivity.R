my_coloc_sensitivity <- function(obj, rule = "H4 > 0.8", trait1_title = "trait 1", trait2_title = "trait 2",
                              dataset1 = NULL, dataset2 = NULL,
                              npoints = 100, suppress_messages = FALSE,
                              row = 1, data_check_trait1 = NULL, data_check_trait2 = NULL) {
  # check ====
  stopifnot("list" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if (rule == "")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)", "PP.\\1.abf", rule, perl = TRUE)
  
  ## massage results object
  results <- obj$results
  ## multiple signals?
  multiple <- FALSE
  if (data.table::is.data.table(obj$summary)) { # we're not in coloc.abf anymore
    if (!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ", nrow(obj$summary))
    pp <- unlist(c(obj$summary[row, grep("PP|nsnp", names(obj$summary)), with = FALSE]))
    if (paste0("SNP.PP.H4.row", row) %in% names(results)) {
      multiple <- TRUE
      results[["SNP.PP.H4"]] <- results[[paste0("SNP.PP.H4.row", row)]]
    }
    if (paste0("z.df1.row", row) %in% names(results)) { # might be passed here or in separate dataset objects
      results[["z.df1"]] <- results[[paste0("z.df1.row", row)]]
      results[["z.df2"]] <- results[[paste0("z.df2.row", row)]]
    } else {
      pp <- unlist(c(obj$summary[row, grep("PP|nsnp", names(obj$summary)), with = FALSE]))
    }
  } else {
    pp <- obj$summary
  }
  ## need to add z score from datasets?
  if (!is.null(dataset1) && !is.null(dataset2)) {
    df1 <- with(dataset1, data.table(snp = snp, position = position, z.df1 = beta / sqrt(varbeta)))
    df2 <- with(dataset2, data.table(snp = snp, position = position, z.df2 = beta / sqrt(varbeta)))
    df <- merge(df1, df2, by = c("snp", "position"), all = TRUE)
    results <- merge(results, df, by = "snp")
  }
  
  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp), eval(parse(text = rule))) }
  pass.init <- check(pp)
  if (!suppress_messages) {
    message("Results ", if (check(pp)) { "pass" } else { "fail" }, " decision rule ", rule.init)
  }
  
  testp12 <- 10^seq(log10(p1 * p2), log10(min(p1, p1)), length.out = npoints)
  testH <- coloc_prior.snp2hyp(pp["nsnps"], p12 = testp12, p1 = p1, p2 = p2)
  testpp <- as.data.frame(coloc_prior.adjust(summ = pp, newp12 = testp12, p1 = p1, p2 = p2, p12 = p12))
  colnames(testpp) <- gsub("(H.)", "PP.\\1.abf", colnames(testpp), perl = TRUE)
  pass <- check(testpp)
  w <- which(pass)
  
  # plots ====
  # locuszoom plots ====
  plot_locuscompare <- locuscomparer_make_plot(data_coloc_exposure = data_check_trait1,
                                               data_coloc_outcome = data_check_trait2,
                                               SNP_causal_exposure = SNP_causal_exposure$snp,
                                               trait1_title = trait1_title, trait2_title = trait2_title)
  
  # Manhattan plots ====
  if ("z.df1" %in% colnames(results) && "z.df2" %in% colnames(results)) {
    plot_manhattan <- cowplot::plot_grid(
      coloc_manh.plot(results, 1),
      coloc_manh.plot(results, 2),
      ncol = 1, nrow = 2)
  } else {
    warning("Please supply dataset1 and dataset2 to view Manhattan plots.")
  }
  
  # Probability plots ====
  m <- list(testH, as.matrix(testpp))
  ti <- list("Prior probabilities", "Posterior probabilities")
  plot_probabilities <- list()
  for (i in 1:2) {
    # Set max y axis value
    ym <- if (i == 1) max(m[[i]][, -1]) else max(m[[i]])
    
    # Data
    df <- data.frame(
      testp12 = testp12,
      m = m[[i]],
      pass = pass
    ) %>%
      tidyr::pivot_longer(
        cols = if (i == 1) c("m.H0", "m.H1", "m.H2", "m.H3", "m.H4")
        else c("m.PP.H0.abf", "m.PP.H1.abf", "m.PP.H2.abf", "m.PP.H3.abf", "m.PP.H4.abf"),
        names_to = "name", values_to = "value"
      ) %>%
      dplyr::mutate(name = forcats::fct_recode(name,
                                               H0 = if (i == 1) "m.H0" else "m.PP.H0.abf",
                                               H1 = if (i == 1) "m.H1" else "m.PP.H1.abf",
                                               H2 = if (i == 1) "m.H2" else "m.PP.H2.abf",
                                               H3 = if (i == 1) "m.H3" else "m.PP.H3.abf",
                                               H4 = if (i == 1) "m.H4" else "m.PP.H4.abf"))
    
    # Check data
    df <- df %>%
      dplyr::filter(!is.na(testp12) & !is.na(value)) %>%
      dplyr::filter(testp12 >= 1e-8 & testp12 <= 1e-4)
    p12 <- ifelse(p12 >= 1e-8 & p12 <= 1e-4, p12, NA)
    
    # Initialize the plot
    plot_probabilities[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = testp12, y = value, colour = name)) +
      ggplot2::scale_x_log10(breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4)) +
      ggplot2::scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::labs(x = "p12", y = "Prob",
                    title = ti[[i]],
                    subtitle = paste("shaded region:", rule.init)) +
      ggplot2::geom_vline(xintercept = p12, linetype = "dashed", color = "gray") +
      ggplot2::annotate("text", x = p12, y = 0.5, label = "results", angle = 90, color = "gray40") +
      cowplot::theme_cowplot() +
      ggplot2::theme(legend.position = "right", legend.title = ggplot2::element_blank())
    
    # Add rectangle if pass condition is met
    if (any(df$pass)) {
      # Calculate x-axis limits for the rectangle
      xleft <- min(df$testp12[df$pass])
      xright <- max(df$testp12[df$pass])
      
      plot_probabilities[[i]] <- plot_probabilities[[i]] +
        ggplot2::geom_rect(ggplot2::aes(xmin = xleft, xmax = xright, ymin = 0, ymax = 1),
                           fill = "#E2FFBF", colour = NA, alpha = 0.005)
    }
    
    # Add the points on top of the rectangle
    plot_probabilities[[i]] <- plot_probabilities[[i]] +
      ggplot2::geom_point()
  }
  # combined
  plot_probabilities <- cowplot::plot_grid(plot_probabilities[[1]],
                                           plot_probabilities[[2]],
                                           ncol = 1, nrow = 2)
  
  # combined plot ====
  plot_sensitivity <- cowplot::plot_grid(
    plot_locuscompare,
    plot_manhattan,
    plot_probabilities,
    nrow = 1, rel_widths = c(1,0.6,0.6)
  )
  
  plot_sensitivity
  
  return(plot_sensitivity)
  invisible(cbind(testpp, p12 = testp12, pass = pass))
}
