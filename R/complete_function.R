#' End-to-End Causal Scan (multi-outcome, robust validation, and binary pairs)
#'
#' @description
#' Runs the full pipeline on \code{df}: (i) column preprocessing (coercion to
#' factor or numeric with rare-level handling), (ii) multi-outcome scan via
#' \code{scan_all_outcomes_complete()}, (iii) robust validation with bootstrap +
#' permutation via \code{robust_scan_all_outcomes()}, and (iv) optional binary
#' pairwise direction matrix via \code{cor_forest_matrix_robust_perm()}.
#' Optionally draws a causal graph that combines robust multi-outcome and binary evidence.
#'
#' @param df A data frame with variables to analyze (predictors + outcomes).
#' @param n_boot Integer. Number of bootstraps (\code{B}) for robust/binary stages. Default \code{100}.
#' @param n_perm Integer. Number of permutations (\code{P}) for robust/binary stages. Default \code{50}.
#' @param alpha Numeric in (0, 1). Significance level for permutation tests. Default \code{0.05}.
#' @param ntree Integer. Number of trees for Random Forest fits. Default \code{500}.
#' @param seed Optional integer random seed.
#' @param n_cores Integer. Parallel cores for \code{parallel::mclapply()}.
#'   Default \code{parallel::detectCores() - 1}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#' @param always_predictors Optional character vector; column names forced to be considered as predictors.
#' @param which One of \code{"all"}, \code{"robust"}, \code{"binary"}:
#'   controls which stages of the pipeline are executed.
#'   If \code{df} has exactly 2 columns, the mode is coerced to \code{"binary"}.
#' @param categorical_thr Numeric. Threshold for categorical variables in \code{evaluate_importance()}.
#'   Default \code{35}.
#' @param quantitative_thr Numeric. Threshold for quantitative variables in \code{evaluate_importance()}.
#'   Default \code{40}.
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"}; strategy used by
#'   \code{evaluate_importance()}.
#' @param prob Numeric. Probability cutoff used when \code{importance_method = "net_clust"}.
#'   Default \code{0.75}.
#' @param plot Logical. If \code{TRUE}, build and print a combined causal graph. Default \code{TRUE}.
#' @param curved See \code{\link{draw_dag}}: which edges to draw as arcs (\code{NULL}, logical vector,
#'   character keys \code{"X->Y"}, or \code{data.frame(from,to)}).
#' @param layout Graph layout name passed to \code{ggraph}. Examples: \code{"auto"}, \code{"kk"},
#'   \code{"fr"}, \code{"sugiyama"}, \code{"linear"}. Default \code{"auto"}.
#' @param pad Numeric padding added around computed x/y ranges for the plot. Default \code{0.4}.
#' @param arrow_len_pt Arrow length (points) for directed edges. Default \code{8}.
#' @param end_cap_mm End cap radius (millimeters) for edge arrows. Default \code{8}.
#' @param linewidth Edge line width. Default \code{0.75}.
#' @param node_size Node point size. Default \code{20}.
#' @param node_stroke Node point stroke width. Default \code{1}.
#' @param strength_curved Curvature strength for curved edges (passed to \code{ggraph::geom_edge_arc2()}).
#'   Non-curved edges use 0. Default \code{0.6}.
#'
#' @details
#' \strong{Preprocessing.} Columns are preprocessed to stabilize the scans:
#' (a) low-cardinality numeric vectors are converted to factors; (b) high-cardinality
#' non-numeric vectors are coerced to numeric when \code{can_coerce_numeric()} returns \code{TRUE}
#' (dropping rows that cannot be safely coerced), otherwise they are factored with rare-level
#' grouping and dominance checks; (c) strongly dominant categorical columns may be binarized
#' or dropped. During this stage, some columns may be appended to \code{always_predictors}.
#'
#' \strong{Scanning.} If \code{which %in% c("all","robust")}, \code{scan_all_outcomes_complete()}
#' runs first, followed by \code{robust_scan_all_outcomes()} which applies bootstrap (\code{B = n_boot})
#' and permutation tests (\code{P = n_perm}) at level \code{alpha}.
#'
#' \strong{Binary pairs.} If \code{which %in% c("all","binary")}, the binary direction matrix is
#' computed via \code{cor_forest_matrix_robust_perm()} with its own bootstrap/permutation routine.
#'
#' \strong{Graph.} When \code{plot = TRUE}, a combined edge set is built by merging robust multi-outcome
#' relations (labeled \code{"complex"}) and significant binary pairs (labeled \code{"pairwise"}).
#' Overlaps are labeled \code{"both"}. The graph is drawn with \code{\link{draw_dag}} using the provided
#' layout/curvature/appearance settings.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{scan_res}: results from \code{scan_all_outcomes_complete()} (or \code{NULL});
#'   \item \code{robust}: results from \code{robust_scan_all_outcomes()} (or \code{NULL});
#'   \item \code{binary}: list of matrices from \code{cor_forest_matrix_robust_perm()} (or \code{NULL});
#'   \item \code{graph}: a \strong{ggplot} object returned by \code{draw_dag()} (or \code{NULL});
#'   \item \code{res_all}: combined edge data.frame used for plotting (from/pairwise/complex/both).
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' out <- complete_function(
#'   df,
#'   n_boot = 50, n_perm = 200, ntree = 400,
#'   which = "all", verbose = TRUE,
#'   plot = TRUE, layout = "kk"
#' )
#'
#' # Combined graph object (ggplot):
#' out$graph
#'
#' # Combined edges used for plotting:
#' out$res_all
#' }
#'
#' @seealso
#' \code{\link{scan_all_outcomes_complete}},
#' \code{\link{robust_scan_all_outcomes}},
#' \code{\link{cor_forest_matrix_robust_perm}},
#' \code{\link{can_coerce_numeric}},
#' \code{\link{draw_dag}}
#'
#' @importFrom parallel detectCores
#' @importFrom stats na.omit
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
                               which = c("all","robust","binary"),
                               categorical_thr = 35,
                               quantitative_thr = 40,
                               importance_method =c("fixed","neg_exp","net_clust"),
                               prob = 0.75,
                               #plot parameters
                               plot = TRUE,
                               curved = NULL,
                               layout = "auto",
                               pad = 0.4,
                               arrow_len_pt = 8,
                               end_cap_mm   = 8,
                               linewidth    = 0.75,
                               node_size    = 20,
                               node_stroke  = 1,
                               strength_curved = 0.6
) {
  # Match mode argument
  which <- match.arg(which)

  # Force binary mode if only 2 variables
  if(ncol(df) == 2){
    which <- "binary"
  if(verbose) cat("\nOnly 2 variables present in the data set. It will run only the binary analysis.\n")
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

  res_rob <- data.frame()
  for(rel in names(robust_scan)){
    if(robust_scan[[rel]]$significant){
      parts <- strsplit(rel, " → ")[[1]]
      preds <- trimws(strsplit(parts[1], "\\+")[[1]])
      out <- trimws(parts[2])
      res_rob <- rbind(res_rob,data.frame(from = preds, to = out, custom_color = 'complex'))
    }
  }

  res_pair <- which(binary_res$significant == 1, arr.ind = TRUE)
  if(!nrow(res_pair) == 0){
    res_pair <- data.frame(
      from = rownames(binary_res$significant)[res_pair[,1]],
      to   = colnames(binary_res$significant)[res_pair[,2]],
      custom_color = 'pairwise')
  }

  res_all <- rbind(res_pair, res_rob)

  if(nrow(res_all) != 0){
    res_all$combo <- paste(res_all$from, res_all$to)

    tab <- table(res_all$combo)

    res_all$custom_color <- ifelse(
      tab[res_all$combo] > 1, "both", res_all$custom_color
    )

    res_all <- res_all[!duplicated(res_all$combo), c("from","to","custom_color")]
  }

  # --- Plot causal graph if requested ---
  g <- NULL
  if (plot & nrow(res_all) != 0) {
    if (is.null(robust_scan) && is.null(binary_res)) {
      if (verbose) cat("\nℹ️  Nothing to plot (no results produced in this mode).\n")
    } else {
      g <- draw_dag(res_all,
                    curved = curved,
                    layout = layout,
                    pad = pad,
                    arrow_len_pt = arrow_len_pt,
                    end_cap_mm   = end_cap_mm,
                    linewidth    = linewidth,
                    node_size    = node_size,
                    node_stroke  = node_stroke,
                    strength_curved = strength_curved)
    }
    print(g)
  }

  if (verbose) cat("\n✅ Done.\n")

  return(list(
    scan_res = scan_results,
    robust   = robust_scan,
    binary   = binary_res,
    graph    = g,
    res_all = res_all
  ))
}

