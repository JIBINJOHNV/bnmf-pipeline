
prep_z_matrix <- function(z_mat, N_mat,
                          corr_cutoff=0.85,
                          rm_traits=NULL,
                          pval_cutoff=NULL) {
  
    # Given a matrix of z-scores (N_variants x M_traits) and vector of median
    # sample sizes per trait:
    # 1) perform final pre-processing steps before bNMF clustering:
    # trait filtering by p-value, trait pruning based on correlation,
    # and z-score scaling based on sample size
    # 2) expand N x M matrix into N x 2M non-negative matrix
    
    # if no user input for pval_cutoff, use bonferroni definition
    if (is.null(pval_cutoff)){
        pval_cutoff <- 0.05 / nrow(z_mat) 
    }
    print(paste0(sum(is.na(z_mat)), " missing values before pval pruning."))
    
    df_traits <- data.frame(trait = colnames(z_mat))
    
    df_traits_filtered <- data.frame(trait=as.character(),
                            result=as.character(),
                            note=as.character())
    
    # remove manually entered traits
    if (!is.null(rm_traits)) {
        print(sprintf("Removing %i manually entered traits!!!", length(rm_traits)))
        z_mat <- z_mat[, !colnames(z_mat) %in% rm_traits]
        
        df_remove_manual <- data.frame(trait=rm_traits,
                                    result="removed (manual)",
                                    note=NA)
        df_traits_filtered <- rbind(df_traits_filtered, df_remove_manual)
    }
    
    # Filter traits by p-value (min. p-value < 0.05/N_variants)
    print("Filtering traits w/ no pvalues below cutoff...")
    minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
    traits_removed <- colnames(z_mat)[minP_vec >= pval_cutoff]
    df_lowP_vec <- data.frame(minPval=minP_vec[minP_vec >= pval_cutoff]) %>%
        rownames_to_column('trait') %>%
        mutate(result="removed (p-value)") %>% 
        mutate(note=paste("min pval=",format(minPval,scientific=T))) %>%
        dplyr::select(trait, result, note)
    
    df_traits_filtered <- rbind(df_traits_filtered, df_lowP_vec)
    cat(paste(sprintf("Removing traits with no variant having p < %.3e",pval_cutoff),
                paste(traits_removed,
                    collapse="\n"),sep=":\n"))
    z_mat <- z_mat[, minP_vec < pval_cutoff]

    # Prune traits by correlation (remove traits with Pearson |r| > 0.85)
    
    cat(sprintf("\n\nPrune traits by correlation (remove traits with Pearson |r| > %.2f)\n",
                corr_cutoff))
    trait_cor_mat <- cor(z_mat, use="pairwise.complete.obs")  # Trait-trait correlation matrix
    write.table(trait_cor_mat,"./trait_cor_mat.txt", sep="\t")
    
    # sort by max(z) instead of min(pval)
    remaining_traits <- names(sort(apply(z_mat, 2, max, na.rm=T),decreasing = T))
    print(paste("Initial number of traits:",length(remaining_traits)))
    
    keep_traits <- c()
    df_corr_removed <- c()
    
    while (length(remaining_traits) > 0) {
        # append to list of traits to keep
        keep_traits <- c(keep_traits, remaining_traits[1])
        if (length(remaining_traits)==1) {
        break
        }
        trait_cor_mat <- trait_cor_mat[remaining_traits, remaining_traits]
        print(dim(trait_cor_mat))
        to_remove <- rownames(trait_cor_mat)[abs(trait_cor_mat[, remaining_traits[1]]) >= corr_cutoff]

        # track results in dataframe
        if (length(to_remove)>1) {
        df_corr_removed_tmp <- data.frame(trait=to_remove[to_remove!=remaining_traits[1]],
                                            result="removed (correlation)",
                                            note=paste("correlated w/", remaining_traits[1]))
        df_traits_filtered <- rbind(df_traits_filtered, df_corr_removed_tmp)
        }
        
        if (length(to_remove) > 1) {
        cat(paste(sprintf("Correlated trait being removed for %s:",remaining_traits[1]),
                    paste(to_remove[!to_remove %in% remaining_traits[1]],collapse="\n"),"\n",sep="\n"))
        }
        
        # remove all correlated traits from remaining list
        remaining_traits <- setdiff(
        remaining_traits, 
        to_remove
        )
    }
    
    df_traits_filtered <- df_traits_filtered %>%
        right_join(df_traits, by="trait") %>%
        mutate(result = ifelse(is.na(result), "trait kept", result))

    pruned_traits <- df_traits_filtered %>%
        filter(result!="trait kept") %>%
        pull(trait)
    cat(paste("Traits removed in pruning process:", 
                paste(pruned_traits, collapse="\n"),sep = "\n"))
    
    cat(paste0("\nNumber of remaining traits: ",length(keep_traits),"\n"))
    
    # cat(paste("Remaing traits:", 
    #           paste(keep_traits, collapse="\n"),sep = "\n")) 
    print(paste0(sum(is.na(z_mat)), " missing values before corr pruning."))
    z_mat <- z_mat[, keep_traits]

    # Adjust z-scores by sample size for each variant-trait combo
    # i.e. (z = z / sqrt(medN) * mean(sqrt(medN_all_traits)))
    cat("\n\n")
    print("Performing sample size adjustment...")
    medN_vec <- apply(N_mat[, colnames(z_mat)], 2, median, na.rm=T)
    z_mat <- z_mat / sqrt(N_mat[, colnames(z_mat)]) * mean(sqrt(medN_vec))
    
    # Replace missing values with zero
    print("Replacing remaining missing values with zero...")
    print(paste0(sum(is.na(z_mat)), " missing values were replaced."))
    z_mat[is.na(z_mat)] <- 0
    
    # Expand into N x 2M non-negative matrix
    print("Expanding z-score matrix into non-negative matrix (N-variants x 2M-traits)...")
    z_mat_pos <- z_mat
    z_mat_pos[z_mat_pos < 0] <- 0
    colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
    z_mat_neg <- -z_mat
    z_mat_neg[z_mat_neg < 0] <- 0
    colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
    final_z_mat <- cbind(z_mat_pos, z_mat_neg)
    
    output <- list(final_z_mat=final_z_mat,
                    df_traits=df_traits_filtered)
}



save_bnmf_results <- function(bnmf_reps,
                              main_gwas_id,
                              output_folder) {

  # -------------------------------
  # 1. Input validation
  # -------------------------------
  if (!is.list(bnmf_reps) || length(bnmf_reps) == 0) {
    stop("bnmf_reps must be a non-empty list (output of run_bNMF).")
  }

  if (!is.character(main_gwas_id) || nchar(main_gwas_id) == 0) {
    stop("main_gwas_id must be a non-empty character string.")
  }

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  message("Saving bNMF results for: ", main_gwas_id)
  message("Output folder: ", normalizePath(output_folder))

  # -------------------------------
  # 2. Run-Level Metrics
  # -------------------------------
  run_metrics <- data.frame(
    run_id = seq_along(bnmf_reps),

    final_active = sapply(bnmf_reps, function(x)
      as.numeric(tail(unlist(x$n.active), 1))),

    iterations = sapply(bnmf_reps, function(x)
      x$iterations),

    converged = sapply(bnmf_reps, function(x)
      x$converged),

    final_evidence = sapply(bnmf_reps, function(x)
      as.numeric(tail(unlist(x$n.evid), 1))),

    final_error = sapply(bnmf_reps, function(x)
      as.numeric(tail(unlist(x$n.error), 1))),

    final_lambda_mean = sapply(bnmf_reps, function(x)
      mean(tail(unlist(x$n.lambda), 1))),

    final_lambda_sd = sapply(bnmf_reps, function(x)
      sd(tail(unlist(x$n.lambda), 1)))
  )

  run_metrics_file <- file.path(output_folder,
    paste0(main_gwas_id, "_bnmf_run_metrics.csv"))

  write.csv(run_metrics, run_metrics_file, row.names = FALSE)

  # -------------------------------
  # 3. K-Level Summary
  # -------------------------------
  k_summary <- run_metrics |>
    dplyr::group_by(final_active) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      proportion = n_runs / nrow(run_metrics),
      mean_evidence = mean(final_evidence),
      sd_evidence = sd(final_evidence),
      mean_error = mean(final_error),
      sd_error = sd(final_error),
      mean_iterations = mean(iterations),
      .groups = "drop"
    ) |>
    dplyr::arrange(desc(proportion))

  k_summary_file <- file.path(output_folder,
    paste0(main_gwas_id, "_bnmf_K_summary.csv"))

  write.csv(k_summary, k_summary_file, row.names = FALSE)

  # -------------------------------
  # 4. Factor Strength Table
  # -------------------------------
  factor_metrics <- do.call(rbind,
    lapply(seq_along(bnmf_reps), function(i) {

      W <- bnmf_reps[[i]]$W
      H <- bnmf_reps[[i]]$H
      k_eff <- run_metrics$final_active[i]

      if (k_eff > 0) {
        data.frame(
          run_id = i,
          factor_id = 1:k_eff,
          W_norm = apply(W[, 1:k_eff, drop = FALSE], 2,
                         function(x) sqrt(sum(x^2))),
          H_norm = apply(H[1:k_eff, , drop = FALSE], 1,
                         function(x) sqrt(sum(x^2)))
        )
      } else {
        NULL
      }
    })
  )

  factor_metrics_file <- file.path(output_folder,
    paste0(main_gwas_id, "_bnmf_factor_metrics.csv"))

  write.csv(factor_metrics, factor_metrics_file, row.names = FALSE)

  # -------------------------------
  # 5. Cluster Loadings Table
  # -------------------------------
  cluster_loadings <- do.call(rbind,
    lapply(seq_along(bnmf_reps), function(i) {

      H <- bnmf_reps[[i]]$H
      k_eff <- run_metrics$final_active[i]

      if (k_eff > 0) {
        df <- as.data.frame(H[1:k_eff, , drop = FALSE])
        df$run_id <- i
        df$factor_id <- 1:k_eff
        df
      } else {
        NULL
      }
    })
  )

  cluster_loadings_file <- file.path(output_folder,
    paste0(main_gwas_id, "_bnmf_cluster_loadings.csv"))

  write.csv(cluster_loadings, cluster_loadings_file, row.names = FALSE)

  message("bNMF result tables saved successfully.")

  # -------------------------------
  # 6. Return objects invisibly
  # -------------------------------
  invisible(list(
    run_metrics = run_metrics,
    k_summary = k_summary,
    factor_metrics = factor_metrics,
    cluster_loadings = cluster_loadings
  ))
}
