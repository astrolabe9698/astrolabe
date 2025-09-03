#source('/mnt/sannaLAB-Temp/mathronno/Astrolabio/final_chapter/Package/scan_utilities.R')

#' Robust Validation of Causal Scan Results (Bootstrap + Permutation)
#'
#' Validates each relation found by a prior scan (e.g. \code{scan_all_outcomes_complete})
#' using (1) bootstrap repetitions on real data to count how often the relation (or
#' a subset with same outcome) reappears, and (2) an empirical permutation test
#' with embedded bootstrap to compute \eqn{p}-values.
#'
#' @param df Data frame with all variables.
#' @param scan_results List of results returned by \code{scan_all_outcomes_complete()}.
#' @param seed Optional integer random seed.
#' @param ntree Integer. Number of trees for Random Forest during re-scans (default \code{500}).
#' @param n_boot Integer. Bootstrap repetitions on real data (default \code{300}).
#' @param n_perm Integer. Number of permutations (each with embedded bootstrap) (default \code{30}).
#' @param alpha Numeric. Significance level for the empirical \eqn{p}-value (default \code{0.05}).
#' @param n_cores Integer. Number of parallel cores for \code{mclapply} (default \code{parallel::detectCores() - 1}).
#' @param verbose Logical. Print progress (default \code{TRUE}).
#' @param always_predictors Optional character vector of predictors to keep available.
#' @param categorical_thr Numeric. Threshold for categorical variables in \code{evaluate_importance()} (default \code{35}).
#' @param quantitative_thr Numeric. Threshold for quantitative variables in \code{evaluate_importance()} (default \code{40}).
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"} for \code{evaluate_importance()}.
#' @param prob Numeric. Probability cutoff for neural network inside \code{evaluate_importance()} when \code{importance_method="net_clust"}.
#'
#' @return A named list (one entry per tested relation) with elements:
#' \itemize{
#'   \item \code{predictors}, \code{outcome};
#'   \item \code{bootstrap_freq}: total bootstrap hits (partial-match counting);
#'   \item \code{max_freq}: theoretical maximum hits (\code{n_boot * length(predictors)});
#'   \item \code{freq_perm}: vector of permutation bootstrap counts (or \code{NA} if skipped);
#'   \item \code{p_empirical}: empirical \eqn{p}-value;
#'   \item \code{significant}: logical flag (\code{p_empirical < alpha}).
#' }
#'
#' @details
#' Bootstrap step: resamples rows with replacement (and shuffles columns order),
#' reruns \code{scan_all_outcomes_complete()} on the subset of variables for the
#' target relation, and counts partial matches (same outcome and predictors included).
#' Permutation step: shuffles all variables within the subset and repeats the
#' bootstrap to build the null distribution of the count; computes \eqn{p}-value.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' val <- robust_scan_all_outcomes(
#'   df, scan_results,
#'   ntree = 500, n_boot = 200, n_perm = 50, alpha = 0.05, n_cores = 4
#' )
#' }
#'
#' @seealso \code{\link{scan_all_outcomes_complete}}, \code{\link{evaluate_importance}}
#' @importFrom parallel mclapply detectCores
#' @importFrom stats quantile
#' @export
robust_scan_all_outcomes <- function(df,
                                     scan_results,
                                     seed = NULL,
                                     ntree = 500,
                                     n_boot = 300,
                                     n_perm = 30,
                                     alpha = 0.05,
                                     n_cores = parallel::detectCores() - 1,
                                     verbose = TRUE,
                                     always_predictors = NULL,
                                     categorical_thr = 35,
                                     quantitative_thr = 40,
                                     importance_method = c("fixed","neg_exp","net_clust"),
                                     prob = 0.75) {

  # ── Helper cross-platform parallel apply ─────────────────────────
  parallel_lapply <- function(X, FUN, n_cores = 1) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      # Assicurarsi che le funzioni siano disponibili ai cluster
      parallel::clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
      parallel::parLapply(cl, X, FUN)
    } else {
      parallel::mclapply(X, FUN, mc.cores = n_cores)
    }
  }

  # ── Set seed ────────────────────────────────
  if(!is.null(seed)) set.seed(seed)
  final_results <- list()

  # ── Helper: count partial matches ──────────
  count_partial_hits <- function(decisions, predictors, outcome) {
    if (length(decisions) == 0 || all(is.na(decisions))) return(0)
    total <- 0
    for (dec in decisions) {
      parts <- strsplit(dec, " → ")[[1]]
      preds <- trimws(strsplit(parts[1], "\\+")[[1]])
      out <- trimws(parts[2])
      if (out != outcome) next
      if (!all(preds %in% predictors) || length(preds) > length(predictors)) return(0)
      total <- total + sum(predictors %in% preds)
    }
    total
  }

  # ── Verbose header ────────────────────────────────
  if (verbose) {
    cat("\n────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
    cat("🧪 Robust validation of scan results (robust_scan_all_outcomes)\n")
    cat(sprintf("   • Relations to test: %d\n", length(scan_results)))
    cat(sprintf("   • Settings: ntree=%d, B=%d (bootstrap), P=%d (permutations), α=%.2f, cores=%d, seed=%s\n",
                ntree, n_boot, n_perm, alpha, n_cores, ifelse(is.null(seed), "NULL", as.character(seed))))
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
  }
  # ── Loop su ogni relazione ─────────────────
  for (i in seq_along(scan_results)) {
    decision <- scan_results[[i]]$drill_down_decision
    if (is.null(decision)) next

    parts <- strsplit(decision, " → ")[[1]]
    predictors <- trimws(strsplit(parts[1], "\\+")[[1]])
    outcome <- trimws(parts[2])
    all_vars <- c(predictors, outcome)
    full_relation <- paste(decision)

    if (verbose) cat(sprintf("\n[%d/%d] Testing relation: %s\n", i, length(scan_results), full_relation))

    # --- 1. Bootstrap sui dati reali ---
    if (verbose) cat(sprintf("   • Step 1 — Bootstrap on real data (B=%d)\n", n_boot))
    boot_rels <- parallel_lapply(1:n_boot, function(b) {
      if(!is.null(seed)) set.seed(seed + b)
      df_boot <- df[sample(1:nrow(df), replace = TRUE), all_vars]
      df_boot <- df_boot[, sample(1:ncol(df_boot))]
      res <- scan_all_outcomes_complete(remove_outliers(df_boot), seed = seed, verbose = FALSE,
                                        ntree = ntree, always_predictors = always_predictors,
                                        categorical_thr = categorical_thr,
                                        quantitative_thr = quantitative_thr,
                                        importance_method = importance_method,
                                        prob = prob)
      if (is.null(res) || all(sapply(res, is.null))) return(NA)
      decs <- character()
      for (j in seq_along(res)) {
        if (is.null(res[[j]])) next
        dec_boot <- if (!is.null(res[[j]]$drill_down_decision)) res[[j]]$drill_down_decision else res[[j]]$outer_layer_decision
        if (!is.null(dec_boot)) decs <- c(decs, dec_boot)
      }
      decs
    }, n_cores = n_cores)

    freq <- sum(sapply(boot_rels, function(decs) count_partial_hits(decs, predictors, outcome)))

    if (verbose) {
      cat(sprintf("   ✅ Bootstrap frequency: %d/%d (max=%d)\n",
                  round(freq/length(predictors)), n_boot, n_boot ))
      cat(sprintf("   • Step 2 — Permutation test (P=%d) with embedded bootstrap (B=%d)\n", n_perm, n_boot))
      cat(sprintf("     Shuffling all variables within subset: %s\n", paste(all_vars, collapse = ", ")))
    }

    # --- 2. Permutation test ---
    if (freq != 0) {
      freq_perm <- numeric(n_perm)
      for (perm in 1:n_perm) {
        if (verbose) { cat(sprintf("\r     Permutation %d/%d", perm, n_perm)); flush.console() }
        df_perm <- as.data.frame(lapply(df[, all_vars], sample))
        perm_boots <- parallel_lapply(1:n_boot, function(b) {
          if(!is.null(seed)) set.seed(seed + perm*1000 + b)
          df_boot <- df_perm[sample(1:nrow(df_perm), replace = TRUE), ]
          df_boot <- df_boot[, sample(1:ncol(df_boot))]
          res <- scan_all_outcomes_complete(remove_outliers(df_boot), seed = seed, verbose = FALSE,
                                            ntree = ntree, always_predictors = always_predictors,
                                            categorical_thr = categorical_thr,
                                            quantitative_thr = quantitative_thr,
                                            importance_method = importance_method,
                                            prob = prob)
          if (is.null(res) || all(sapply(res, is.null))) return(NA)
          decs <- character()
          for (j in seq_along(res)) {
            if (is.null(res[[j]])) next
            dec_boot <- if (!is.null(res[[j]]$drill_down_decision)) res[[j]]$drill_down_decision else res[[j]]$outer_layer_decision
            if (!is.null(dec_boot)) decs <- c(decs, dec_boot)
          }
          decs
        }, n_cores = n_cores)
        freq_perm[perm] <- sum(sapply(perm_boots, function(decs) count_partial_hits(decs, predictors, outcome)))
      }
      if (verbose) cat("\n")

      # Compute empirical p-value
      p_empirical <- (sum(freq_perm >= freq) + 0.1) / (n_perm + 0.1)

      if (verbose) {
        sig_flag <- ifelse(p_empirical < alpha, "SIGNIFICANT ✳", "n.s.")
        q <- quantile(freq_perm/length(predictors), probs = c(0.25, 0.5, 0.75))
        cat(sprintf("   📊 Result: Distribution of permuted bootstrap (B=%d) Q1=%.0f | Q2=%.0f | Q3=%.0f | p_empirical=%.4g | %s (α=%.2f)\n",
                    n_boot,q[1], q[2], q[3], p_empirical, sig_flag, alpha))
        cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
      }

    } else {
      # Skip permutations if bootstrap freq is zero
      p_empirical <- 1
      freq_perm <- NA
      if (verbose) {
        cat("\nSkipping permutations since bootstrap frequency is 0.\n")
      }
    }

    is_significant <- p_empirical < alpha
    final_results[[full_relation]] <- list(
      predictors = predictors,
      outcome = outcome,
      p_empirical = p_empirical,
      bootstrap_freq = freq,
      max_freq = n_boot * length(predictors),
      freq_perm = freq_perm,
      significant = is_significant
    )
  }

  # ── Verbose summary ────────────────────────────────
  if (verbose) {
    n_sig <- sum(vapply(final_results, function(x) isTRUE(x$significant), logical(1)))
    cat(sprintf("\n✅ Robust scan complete. Significant relations: %d\n", n_sig))
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
  }


  return(final_results)
}



#' Robust Pairwise Direction Matrix via Bootstrap + Permutation
#'
#' For every unordered variable pair \eqn{(X, Y)}, infers the preferred direction
#' (\eqn{X \to Y} or \eqn{Y \to X}) on real data, then estimates robustness by
#' bootstrap frequency and an empirical permutation test. Aggregates results into
#' four matrices.
#'
#' @param df Data frame with all variables.
#' @param B Integer. Number of bootstrap repetitions on real data (default \code{30}).
#' @param P Integer. Number of permutations (each with embedded bootstrap) (default \code{30}).
#' @param seed Optional integer random seed.
#' @param ntree Integer. Trees for Random Forest in inner scans (default \code{300}).
#' @param n_cores Integer. Parallel cores for \code{mclapply} (default \code{parallel::detectCores() - 1}).
#' @param alpha Numeric. Significance level for empirical \eqn{p}-values (default \code{0.05}).
#' @param verbose Logical. Print progress (default \code{TRUE}).
#' @param always_predictors Optional character vector of predictors to keep available (pairs where both are in this set are skipped).
#' @param categorical_thr,quantitative_thr Thresholds for \code{evaluate_importance()}.
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"}.
#' @param prob Numeric. Probability cutoff for neural network when \code{importance_method="net_clust"}.
#'
#' @return A list of four matrices (dimension = \code{p × p}, row/col names = variable names):
#' \itemize{
#'   \item \code{real}: binary adjacency (1 if direction selected on real data);
#'   \item \code{freq}: bootstrap hit counts for the selected direction;
#'   \item \code{significant}: binary adjacency after permutation test (\code{p < alpha});
#'   \item \code{pval}: matrix of empirical \eqn{p}-values (rounded).
#' }
#'
#' @details
#' For each pair, runs \code{scan_all_outcomes_complete()} on the 2D subset and
#' keeps only truly binary decisions. Bootstrap counts how often the same direction
#' reappears; permutations shuffle the target outcome to form a null distribution.
#'
#' @examples
#' \dontrun{
#' mats <- cor_forest_matrix_robust_perm(
#'   df, B = 50, P = 100, ntree = 300, n_cores = 4, alpha = 0.05
#' )
#' image(mats$significant)  # quick look at significant directions
#' }
#'
#' @seealso \code{\link{scan_all_outcomes_complete}}, \code{\link{robust_scan_all_outcomes}}
#' @importFrom parallel mclapply detectCores
#' @importFrom stats quantile
#' @export
cor_forest_matrix_robust_perm <- function(df,
                                          B = 30,
                                          P = 30,
                                          seed = NULL,
                                          ntree = 300,
                                          n_cores = parallel::detectCores() - 1,
                                          alpha = 0.05,
                                          verbose = TRUE,
                                          always_predictors = NULL,
                                          categorical_thr = 35,
                                          quantitative_thr = 40,
                                          importance_method = c("fixed","neg_exp","net_clust"),
                                          prob = 0.75) {

  # ── Cross-platform parallel apply ─────────────
  parallel_lapply <- function(X, FUN, n_cores = 1) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      # esporta tutto l'ambiente ai cluster
      parallel::clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
      parallel::parLapply(cl, X, FUN)
    } else {
      parallel::mclapply(X, FUN, mc.cores = n_cores)
    }
  }

  # ── Set seed ─────────────
  if(!is.null(seed)) set.seed(seed)
  vars <- colnames(df)
  p <- length(vars)

  # ── Inizializza matrici ─────────────
  mat_real <- matrix(0, nrow = p, ncol = p, dimnames = list(vars, vars))
  mat_freq <- matrix(0, nrow = p, ncol = p, dimnames = list(vars, vars))
  mat_sig  <- matrix(0, nrow = p, ncol = p, dimnames = list(vars, vars))
  mat_pval <- matrix(NA, nrow = p, ncol = p, dimnames = list(vars, vars))

  combs <- combn(vars, 2, simplify = FALSE)

  if (verbose) {
    cat("\n────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
    cat("🧭 Beginning robust pairwise scan (cor_forest_matrix_robust_perm)\n")
    cat(sprintf("   • Variables: %d | Pairs: %d\n", length(vars), length(combs)))
    cat(sprintf("   • Settings: B=%d (bootstrap), P=%d (permutations), ntree=%d, α=%.2f, cores=%d, seed=%s\n",
                B, P, ntree, alpha, n_cores, ifelse(is.null(seed), "NULL", as.character(seed))))
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
  }

  results <- lapply(combs, function(pair) {
    var1 <- pair[1]; var2 <- pair[2]
    if (verbose) message(sprintf("\n➡️  Pair: %s ↔ %s\n", var1, var2))
    if(var1 %in% always_predictors & var2 %in% always_predictors) {
      if (verbose) cat("   ⏭️  Skipping: both in always_predictors.\n")
      return(NULL)
    }

    df_real <- df[, c(var1, var2)]
    res <- scan_all_outcomes_complete(df_real, ntree = ntree, seed = seed, verbose = FALSE,
                                      always_predictors = always_predictors,
                                      categorical_thr = categorical_thr,
                                      quantitative_thr = quantitative_thr,
                                      importance_method = importance_method,
                                      prob = prob)
    if (length(res) == 0)  {
      if (verbose) cat("   ⚠️  No decision found on real data for this pair — skipping.\n")
      return(NULL)
    }

    dec <- if (!is.null(res[[1]]$drill_down_decision)) res[[1]]$drill_down_decision else res[[1]]$outer_layer_decision
    parts <- strsplit(dec, " → ")[[1]]
    predictors <- trimws(strsplit(parts[1], "\\+")[[1]])
    outcome <- trimws(parts[2])

    if (length(predictors) != 1 || !(outcome %in% c(var1, var2))) {
      if (verbose) cat(sprintf("   ℹ️  Non-binary decision for this pair (%s) — skipping.\n", dec))
      return(NULL)
    }
    from <- predictors
    to <- outcome

    if (verbose) {
      cat(sprintf("   ✔️  Candidate direction: %s → %s\n", from, to))
      cat(sprintf("   🔁 Bootstrapping real data (B=%d)...\n", B))
    }

    freq_boot <- unlist(parallel_lapply(1:B, function(b) {
      if(!is.null(seed)) set.seed(seed + b)
      df_boot <- df[sample(nrow(df), replace = TRUE), c(var1, var2)]
      res_boot <- scan_all_outcomes_complete(remove_outliers(df_boot), seed = seed, ntree = ntree, verbose = FALSE,
                                             always_predictors = always_predictors,
                                             categorical_thr = categorical_thr,
                                             quantitative_thr = quantitative_thr,
                                             importance_method = importance_method,
                                             prob = prob)
      if (length(res_boot) == 0) return(0L)
      dec_boot <- if (!is.null(res_boot[[1]]$drill_down_decision)) res_boot[[1]]$drill_down_decision else res_boot[[1]]$outer_layer_decision
      parts_boot <- strsplit(dec_boot, " → ")[[1]]
      pred_boot <- trimws(strsplit(parts_boot[1], "\\+")[[1]])
      outc_boot <- trimws(parts_boot[2])
      if (identical(pred_boot, from) && outc_boot == to) 1L else 0L
    }, n_cores = n_cores))

    count_from_to <- sum(freq_boot)
    if (verbose) cat(sprintf("   ✅ Bootstrap frequency: %d/%d\n", count_from_to, B))

    # Skip permutation test if frequency too low
    if (count_from_to < 3) {
      if (verbose) cat("   ⚠️  Frequency too low (<3); skipping permutation test for this pair.\n")
      return(list(real = c(from, to),
                  freq = c(from, to, count_from_to),
                  sig = NULL))
    }

    if (verbose) {
      cat(sprintf("   🎲 Permutation test (P=%d): shuffling outcome '%s'...\n", P, to))
    }
    # --- Parallelized permutation test ---
    freq_perm <- unlist(parallel_lapply(1:B, function(b) {
      df_perm <- df
      df_perm[[to]] <- sample(df_perm[[to]])  # shuffle outcome
      df_perm <- df_perm[, c(var1, var2)]

      res_perm <- scan_all_outcomes_complete(remove_outliers(df_perm), seed = seed, ntree = ntree, verbose = FALSE,
                                             always_predictors = always_predictors,
                                             categorical_thr = categorical_thr,
                                             quantitative_thr = quantitative_thr,
                                             importance_method = importance_method,
                                             prob = prob)
      if (length(res_perm) == 0) return(0L)

      # Extract decision for permuted sample
      dec_perm <- if (!is.null(res_perm[[1]]$drill_down_decision)) res_perm[[1]]$drill_down_decision else res_perm[[1]]$outer_layer_decision
      parts_perm <- strsplit(dec_perm, " → ")[[1]]
      pred_perm <- trimws(strsplit(parts_perm[1], "\\+")[[1]])
      outc_perm <- trimws(parts_perm[2])

      if (identical(pred_perm, from) && outc_perm == to) {
        return(1L)
      } else {
        return(0L)
      }
    }, n_cores = n_cores))

    # Empirical p-value
    p_emp <- (sum(freq_perm >= count_from_to) + 0.1) / (P + 0.1)

    if (verbose) {
      sig_flag <- ifelse(p_emp < alpha, "SIGNIFICANT ✳", "n.s.")
      q <- quantile(freq_perm, probs = c(0.25, 0.5, 0.75))

      cat(sprintf("   📊 Result: Distribution of permuted bootstrap (B=%d) Q1=%.0f | Q2=%.0f | Q3=%.0f | p_empirical=%.4g | %s (α=%.2f)\n",
                  B, q[1], q[2], q[3], p_emp, sig_flag, alpha))
    }

    sig <- if (p_emp < alpha) c(from, to) else NULL
    return(list(real = c(from, to), freq = c(from, to, count_from_to), sig = sig, pval = c(from, to, p_emp)))
  })

  # ── Aggregate results ─────────────
  if (verbose) {
    cat("\n🧩 Aggregating matrices...\n")
  }

  for (res in results) {
    if (length(res) == 0) next
    from <- res$real[1]; to <- res$real[2]
    mat_real[from, to] <- 1
    mat_freq[res$freq[1], res$freq[2]] <- as.numeric(res$freq[3])
    if (!is.null(res$sig)) mat_sig[res$sig[1], res$sig[2]] <- 1
    if (!is.null(res$pval)) mat_pval[res$pval[1], res$pval[2]] <- round(as.numeric(res$pval[3]), 2)
  }

  # Verbose summary
  if (verbose) {
    n_sig <- sum(mat_sig)
    cat(sprintf("✅ Done. Significant directed pairs: %d\n", n_sig))
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
  }


  return(list(real = mat_real, freq = mat_freq, significant = mat_sig, pval = mat_pval))
}

