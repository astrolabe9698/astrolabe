#' End-to-End Causal Scan (multi-outcome, robust validation, and binary pairs)
#'
#' Runs the full pipeline on \code{df}: (i) preprocessing of columns (coercion to
#' factor or numeric with rare-level handling), (ii) multi-outcome scan via
#' \code{scan_all_outcomes_complete()}, (iii) robust validation with bootstrap +
#' permutation via \code{robust_scan_all_outcomes()}, and (iv) optional binary
#' pairwise direction matrix via \code{cor_forest_matrix_robust_perm()}.
#' Optionally draws a causal graph combining robust multi-outcome and binary evidence.
#'
#' @param df Data frame with variables to analyze (predictors + outcomes).
#' @param n_boot Integer. Number of bootstraps (\code{B}) for robustness/binary stages (default \code{100}).
#' @param n_perm Integer. Number of permutations (\code{P}) for robustness/binary stages (default \code{50}).
#' @param alpha Numeric in (0,1). Significance level for permutation tests (default \code{0.05}).
#' @param ntree Integer. Number of trees for Random Forest fits (default \code{500}).
#' @param seed Optional integer random seed.
#' @param n_cores Integer. Parallel cores for \code{parallel::mclapply()} (default \code{parallel::detectCores()-1}).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#' @param always_predictors Optional character vector; columns forced to be considered as predictors.
#' @param plot Logical. If \code{TRUE}, build a combined causal graph (default \code{TRUE}).
#' @param which One of \code{"all"}, \code{"robust"}, \code{"binary"}:
#'   controls which parts of the pipeline are executed.
#'   With 2 columns in \code{df} it coerces to \code{"binary"}.
#' @param categorical_thr Numeric. Threshold for categorical variables in \code{evaluate_importance()} (default \code{35}).
#' @param quantitative_thr Numeric. Threshold for quantitative variables in \code{evaluate_importance()} (default \code{40}).
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"}; strategy for \code{evaluate_importance()}.
#' @param prob Numeric. Probability cutoff for the neural network when \code{importance_method="net_clust"} (default \code{0.75}).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{scan_res}: results from \code{scan_all_outcomes_complete()} (or \code{NULL});
#'   \item \code{robust}: results from \code{robust_scan_all_outcomes()} (or \code{NULL});
#'   \item \code{binary}: list of matrices from \code{cor_forest_matrix_robust_perm()} (or \code{NULL});
#'   \item \code{graph}: igraph (or similar) object returned by \code{plot_causal_graph_igraph()} (or \code{NULL});
#'   \item \code{plot}: a recorded plot (from \code{grDevices::recordPlot()}) if \code{plot=TRUE}, else \code{NULL}.
#' }
#'
#' @details
#' Columns are preprocessed as follows: low-cardinality numeric are converted to factors;
#' high-cardinality non-numeric are coerced to numeric when safe (see \code{can_coerce_numeric()}),
#' otherwise factored with rare-level grouping and dominance checks. Some columns may be
#' appended to \code{always_predictors} during this stage to stabilize scanning.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' out <- complete_function(
#'   df, n_boot = 50, n_perm = 200, ntree = 400,
#'   which = "all", verbose = TRUE
#' )
#' out$graph  # combined causal graph object
#' }
#'
#' @seealso
#' \code{\link{scan_all_outcomes_complete}},
#' \code{\link{robust_scan_all_outcomes}},
#' \code{\link{cor_forest_matrix_robust_perm}},
#' \code{\link{can_coerce_numeric}},
#' \code{\link{plot_causal_graph_igraph}}
#'
#' @importFrom parallel detectCores
#' @importFrom stats na.omit
#' @importFrom grDevices recordPlot
#' @export
complete_function <- function(df,
                              n_boot = 100,
                              n_perm = 50,
                              alpha = 0.05,
                              ntree = 500,
                              seed = NULL,
                              n_cores = parallel::detectCores() - 1,
                              verbose = TRUE,
                              always_predictors = NULL,
                              plot = TRUE,
                              which = c("all","robust","binary"),
                              categorical_thr = 35,
                              quantitative_thr = 40,
                              importance_method =c("fixed","neg_exp","net_clust"),
                              prob = 0.75) {

  # Match mode argument
  which <- match.arg(which)

  # Force binary mode if only 2 variables
  if(ncol(df) == 2){
    which <- "binary"
    cat("\nOnly 2 variables present in the data set. It will run only the binary analysis.\n")
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Remove rows with NAs
  if(any(is.na(df))){
    init_row <- nrow(df)
    df <- na.omit(df)
    final_row <- nrow(df)
    cat('\nRemoving', init_row - final_row, 'NAs...\n')
  }

  n <- nrow(df)
  min_per_level <- max(5, floor(0.01 * n))
  rare_thresh <- min_per_level
  dominance_thresh <- 0.95 # dominance threshold

  # --- Preprocess columns: factors, numerics, high-cardinality ---
  for (col in colnames(df)) {
    x <- df[[col]]
    k <- length(unique(x))

    if (!is.numeric(x)) {
      # Factor handling
      if (k <= min(20, ceiling(0.05 * n)) && all(table(x) >= min_per_level)) {
        cat('\nColumn', toupper(col), '→ FACTOR (k =',paste0(k, ')'),'\n')
        df[[col]] <- as.factor(x)
        always_predictors <- c(always_predictors, col)
      } else {
        # Try coercion to numeric
        if (can_coerce_numeric(x)) {
          cat('\nColumn', toupper(col), 'high-cardinality NUMERIC-LOOKING → NUMERIC\n')
          suppressWarnings(df[[col]] <- as.numeric(as.character(x)))
          if(any(is.na(df[[col]]))){
            init_row <- nrow(df)
            df <- na.omit(df)
            final_row <- nrow(df)
            cat('\n   → Removed', init_row - final_row,' rows since they cannot be coerced as numeric\n')
          }
        } else {
          # Factor with rare levels grouped
          cat('\nColumn', toupper(col), 'high-cardinality NON-NUMERIC → FACTOR with rare levels grouped\n')
          x <- as.character(x)
          freq <- table(x)
          rare_levels <- names(freq)[freq < rare_thresh]
          if (length(rare_levels)) x[x %in% rare_levels] <- "OTHER"
          freq_after <- table(x)

          # Dominance check
          if (max(freq_after) / n > dominance_thresh) {
            top_level <- names(which.max(freq_after))
            if (length(freq_after) == 2) {
              cat('  → Binarizing:', top_level, 'vs others\n')
              x <- factor(x == top_level, levels = c(FALSE, TRUE))
              always_predictors <- c(always_predictors, col)
            } else {
              cat('  → Dropping column due to dominance of', top_level, '\n')
              df[[col]] <- NULL
              next
            }
          }
          df[[col]] <- as.factor(x)
          always_predictors <- c(always_predictors, col)
        }
      }
    } else {
      # Numeric columns with few unique values converted to factor
      if (k <= 8 || k <= 0.02 * n) {
        cat('\nColumn', toupper(col), 'numeric with few uniques (k =', paste0(k, ')'), '→ FACTOR\n')
        df[[col]] <- as.factor(x)
        always_predictors <- c(always_predictors, col)
      } else {
        if (verbose) cat('\nColumn', toupper(col), '→ NUMERIC (k =', paste0(k, ')'),'\n')
      }
    }
  }

  # Initialize results
  scan_results <- NULL
  robust_scan  <- NULL
  binary_res   <- NULL

  # Verbose header
  if (verbose) {
    cat("\n────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
    cat("🚦 complete_function\n")
    cat(sprintf("   • Mode: %s\n", which))
    cat(sprintf("   • Settings: ntree=%d, B=%d, P=%d, α=%.2f, seed=%s, cores=%d\n",
                ntree, n_boot, n_perm, alpha,
                ifelse(is.null(seed), "NULL", as.character(seed)),
                n_cores))
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
  }

  # --- SCAN + ROBUST ---
  if (which %in% c("all","robust")) {
    if (verbose) cat("\n🔎 Running initial scan...\n")
    scan_results <- scan_all_outcomes_complete(
      df,
      ntree = ntree,
      verbose = verbose,
      always_predictors = always_predictors,
      seed = seed,
      categorical_thr = categorical_thr,
      quantitative_thr = quantitative_thr,
      importance_method = importance_method,
      prob = prob
    )

    if (verbose) cat("\n🧪 Starting robust scan...\n")
    robust_scan <- robust_scan_all_outcomes(
      df,
      seed = seed,
      scan_results,
      ntree = ntree,
      n_boot = n_boot,
      n_perm = n_perm,
      alpha = alpha,
      n_cores = n_cores,
      verbose = verbose,
      always_predictors = always_predictors,
      categorical_thr = categorical_thr,
      quantitative_thr = quantitative_thr,
      importance_method = importance_method,
      prob = prob
    )
  }

  # --- BINARY scan ---
  if (which %in% c("all","binary")) {
    if (verbose) cat("\n🧰 Starting binary scan...\n")
    binary_res <- cor_forest_matrix_robust_perm(
      df,
      seed = seed,
      B = n_boot,
      P = n_perm,
      verbose = verbose,
      n_cores = n_cores,
      ntree = ntree,
      always_predictors = always_predictors,
      categorical_thr = categorical_thr,
      quantitative_thr = quantitative_thr,
      importance_method = "fixed",
      prob = prob
    )
  }

  # --- Plot causal graph if requested ---
  p <- NULL
  g <- NULL
  if (plot) {
    if (is.null(robust_scan) && is.null(binary_res)) {
      if (verbose) cat("\nℹ️  Nothing to plot (no results produced in this mode).\n")
    } else {
      g <- plot_causal_graph_igraph(robust_scan, binary_res)
      p <- recordPlot()
    }
  }

  if (verbose) cat("\n✅ Done.\n")

  return(list(
    scan_res = scan_results,
    robust   = robust_scan,
    binary   = binary_res,
    graph    = g,
    plot     = p
  ))
}
