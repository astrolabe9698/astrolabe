#' Mixed PCA and Top-Contributing Variables
#'
#' Runs a mixed PCA (\code{PCAmix}) on numeric + categorical data, selects the
#' number of components needed to reach a cumulative explained variance threshold,
#' and extracts the top-\code{top_k} contributing variables for each retained component.
#'
#' @param df Data frame containing numeric and/or categorical variables.
#'   Character columns are coerced to factors.
#' @param ncomp Integer. Maximum number of components to compute (default \code{5}).
#' @param top_k Integer. Number of top contributing variables to keep per component (default \code{5}).
#' @param var_threshold Numeric in (0,1]. Target cumulative explained variance (default \code{0.90}).
#' @param always_predictors Optional character vector (kept for API symmetry; not used internally).
#' @param verbose Logical. If \code{TRUE}, prints progress (default \code{FALSE}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{components}: matrix/data frame of individual component scores;
#'   \item \code{top_vars}: data frame with columns \code{component}, \code{variable}, \code{contribution_percent};
#'   \item \code{contributions}: matrix of variable contributions per component;
#'   \item \code{fit}: the fitted \code{PCAmix} object.
#' }
#'
#' @examples
#' \dontrun{
#' library(PCAmixdata)
#' set.seed(1)
#' df <- data.frame(
#'   x1 = rnorm(200),
#'   x2 = rnorm(200),
#'   g1 = sample(letters[1:3], 200, TRUE),
#'   g2 = sample(c("A","B","C","D"), 200, TRUE)
#' )
#' res <- mixed_pca_topvars(df, ncomp = 4, top_k = 3, var_threshold = 0.8)
#' head(res$top_vars)
#' }
#'
#' @importFrom PCAmixdata PCAmix
#' @export
mixed_pca_topvars <- function(df, ncomp = 5, top_k = 5, var_threshold = 0.90, always_predictors = NULL, verbose = FALSE) {
  stopifnot(is.data.frame(df))  # ensure input is a data frame

  # Identify numeric and categorical columns
  is_num <- vapply(df, is.numeric, logical(1))
  is_cat <- vapply(df, function(x) is.factor(x) || is.character(x), logical(1))

  # Convert character columns to factor if needed
  if (any(is_cat)) df[is_cat] <- lapply(df[is_cat], function(x) if (is.factor(x)) x else factor(x))

  X.quanti <- if (any(is_num)) as.data.frame(df[, is_num, drop = FALSE]) else NULL
  X.quali  <- if (any(is_cat)) as.data.frame(df[, is_cat, drop = FALSE]) else NULL

  # Perform mixed PCA
  fit_all <- PCAmix(X.quanti = X.quanti, X.quali = X.quali, ndim = min(ncomp, ncol(df)), graph = FALSE, rename.level = TRUE)

  # Calculate explained variance and cumulative variance
  var_explained <- fit_all$eig[, 2] / 100
  cum_var <- cumsum(var_explained)

  # Determine how many components reach the variance threshold
  needed_comp <- which(cum_var >= var_threshold)[1]
  if (is.na(needed_comp)) needed_comp <- length(cum_var)
  final_ncomp <- min(needed_comp, ncomp)

  message(sprintf("%d components selected (%.1f%% of explained variance)", final_ncomp, cum_var[final_ncomp] * 100))

  # Refit PCA using only the selected number of components
  fit <- PCAmix(X.quanti = X.quanti, X.quali = X.quali, ndim = final_ncomp, graph = FALSE, rename.level = TRUE)
  comp_scores <- fit$ind$coord  # component scores for individuals

  # Extract contributions for numeric and categorical variables
  contrib_list <- list()
  if (!is.null(X.quanti)) {
    cq <- fit$quanti$contrib
    rownames(cq) <- colnames(X.quanti)
    contrib_list[["num"]] <- cq
  }
  if (!is.null(X.quali)) {
    cc <- fit$quali$contrib
    rownames(cc) <- colnames(X.quali)
    contrib_list[["cat"]] <- cc
  }

  # Merge contributions into a single matrix
  contrib_mat <- if (length(contrib_list) == 2) rbind(contrib_list[[1]], contrib_list[[2]]) else contrib_list[[1]]
  contrib_mat <- contrib_mat[intersect(colnames(df), rownames(contrib_mat)), , drop = FALSE]

  # Identify top contributing variables for each component
  top_list <- lapply(seq_len(final_ncomp), function(j) {
    cj <- contrib_mat[, j]
    ord <- order(cj, decreasing = TRUE)
    take <- head(ord, min(top_k, length(ord)))
    data.frame(
      component = paste0("Dim", j),
      variable  = rownames(contrib_mat)[take],
      contribution_percent = as.numeric(cj[take]),
      stringsAsFactors = FALSE
    )
  })
  top_df <- do.call(rbind, top_list)

  # Return results
  list(
    components = comp_scores,
    top_vars = top_df,
    contributions = contrib_mat,
    fit = fit
  )
}


#' PCA-Guided Causal Scan and Augmentation
#'
#' Preprocesses variables, runs mixed PCA to find top-contributing variables,
#' then re-runs the full causal pipeline (\code{complete_function}) restricted
#' to those variables (optionally re-adding \code{fixed_variables}). Optionally
#' renders a causal graph from the restricted scan.
#'
#' @param df Data frame with variables to analyze.
#' @param ncomp Integer. Maximum number of PCA components (passed to \code{mixed_pca_topvars}).
#' @param top_k Integer. Top variables per component to retain.
#' @param alpha Numeric. Significance level for permutation tests (default \code{0.05}).
#' @param n_boot Integer. Number of bootstraps (default \code{30}).
#' @param n_perm Integer. Number of permutations (default \code{100}).
#' @param ntree Integer. Trees for Random Forest (default \code{500}).
#' @param seed Optional integer seed.
#' @param n_cores Integer. Parallel cores (default \code{parallel::detectCores()-1}).
#' @param which One of \code{"all"}, \code{"robust"}, \code{"binary"}; passed to \code{complete_function}.
#' @param always_predictors Optional character vector of predictors to always include.
#' @param verbose Logical. Verbose output (default \code{TRUE}).
#' @param plot Logical. If \code{TRUE}, build a causal graph from the restricted scan.
#' @param categorical_thr,quantitative_thr Numeric thresholds for \code{evaluate_importance()}.
#' @param fixed_variables Optional character vector of columns to exclude from PCA but re-attach for scanning.
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"} for \code{evaluate_importance()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pca}: result from \code{mixed_pca_topvars};
#'   \item \code{top_variables}: unique variables selected from PCA;
#'   \item \code{scan_topvars}: list returned by \code{complete_function()} on the restricted set;
#'   \item \code{plot}: recorded plot (if \code{plot=TRUE}), else \code{NULL}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(101)
#' res <- pca_scan_and_augment(
#'   df,
#'   ncomp = 5, top_k = 5, which = "all",
#'   n_boot = 30, n_perm = 100, ntree = 500
#' )
#' res$top_variables
#' }
#'
#' @seealso \code{\link{mixed_pca_topvars}}, \code{\link{complete_function}}
#' @importFrom PCAmixdata PCAmix
#' @importFrom grDevices recordPlot
#' @importFrom parallel detectCores
#' @export
pca_scan_and_augment <- function(df,
                                 ncomp = 5,
                                 top_k = 5,
                                 alpha = 0.05,
                                 n_boot = 30, n_perm = 100, ntree = 500,
                                 seed = 101, n_cores = parallel::detectCores()-1,
                                 which = c("all","robust","binary"),
                                 always_predictors = NULL,
                                 verbose = TRUE,
                                 plot = TRUE,
                                 categorical_thr = 35,
                                 quantitative_thr = 40,
                                 fixed_variables = NULL,
                                 importance_method = c("fixed","neg_exp","net_clust")) {
  which <- match.arg(which)

  # Remove rows with missing values
  if(any(as.logical(sum(is.na(df))))){
    init_row <- nrow(df)
    df <- na.omit(df)
    final_row <- nrow(df)
    cat('\nRemoved', init_row - final_row,'rows containing NAs.\n')
  }

  # Preprocess each column to classify as numeric or factor and handle special cases
  for (col in colnames(df)) {
    x <- df[[col]]
    k <- length(unique(x))

    # Drop columns with only one unique value
    if (k == 1) {
      df[[col]] <- NULL
      next
    }

    n <- nrow(df)
    min_per_level <- max(5, floor(0.01 * n))   # Minimum count per factor level
    rare_thresh <- min_per_level                # Threshold for grouping rare levels
    dominance_thresh <- 0.95                    # Threshold for dominant levels

    if (!is.numeric(x)) {
      # Low-cardinality categorical variable → convert to factor
      if (k <= min(20, ceiling(0.05 * n)) && all(table(x) >= min_per_level)) {
        cat('\nColumn', toupper(col), '→ FACTOR (k =', paste0(k, ')'), '\n')
        df[[col]] <- as.factor(x)
        always_predictors <- c(always_predictors, col)
      } else {
        # Try coercing high-cardinality variable to numeric if possible
        if (can_coerce_numeric(x)) {
          cat('\nColumn', toupper(col), 'high-cardinality NUMERIC-LOOKING → NUMERIC\n')
          suppressWarnings(df[[col]] <- as.numeric(as.character(x)))
          if(any(as.logical(sum(is.na(df))))){
            init_row <- nrow(df)
            df <- na.omit(df)
            final_row <- nrow(df)
            cat('\n   → Removed', init_row - final_row, 'rows that could not be coerced to numeric\n')
          }
        } else {
          # High-cardinality non-numeric → factor with rare levels grouped as "OTHER"
          cat('\nColumn', toupper(col), 'high-cardinality NON-NUMERIC → FACTOR with rare levels grouped\n')
          x <- as.character(x)
          freq <- table(x)
          rare_levels <- names(freq)[freq < rare_thresh]
          if (length(rare_levels)) {
            x[x %in% rare_levels] <- "OTHER"
          }
          freq_after <- table(x)

          # Drop column if a single level dominates
          if (max(freq_after) / n > dominance_thresh) {
            top_level <- names(which.max(freq_after))
            if (length(freq_after) == 2) {
              # Binary case → convert to TRUE/FALSE
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
      # Numeric variable
      if (k <= 8 || k <= 0.02 * n) {
        # Few unique numeric values → treat as factor
        cat('\nColumn', toupper(col), 'numeric with few uniques (k =', paste0(k, ')'), '→ FACTOR\n')
        df[[col]] <- as.factor(x)
        always_predictors <- c(always_predictors, col)
      } else {
        if (verbose) cat('\nColumn', toupper(col), '→ NUMERIC (k =', paste0(k, ')'), '\n')
      }
    }
  }

  # Remove fixed variables if specified (used for PCA)
  df_reduced <- df
  if (!is.null(fixed_variables)){
    for (col in fixed_variables){
      df_reduced[[col]] <- NULL
    }
  }

  # Perform mixed PCA and extract top variables contributing to components
  pca <- mixed_pca_topvars(df_reduced, ncomp = ncomp, top_k = top_k)
  comp_df <- as.data.frame(pca$components)
  colnames(comp_df) <- paste0("Dim", seq_len(ncol(comp_df)))

  # Identify variables involved in top components
  tv_involved <- pca$top_vars
  vars_subset <- unique(tv_involved$variable)

  # Re-run the causal scan only on the most important variables
  if (verbose) cat(sprintf("\n🧪 Re-scan on %d important variables...\n", length(vars_subset)))
  df_vars <- df[, intersect(vars_subset, colnames(df)), drop = FALSE]

  # Add fixed variables back for the scan if needed
  if(!is.null(fixed_variables)){
    df_vars <- as.data.frame(cbind(df_vars, df[, fixed_variables]))
    colnames(df_vars)[(length(vars_subset)+1):ncol(df_vars)] <- fixed_variables
  }

  # Run complete causal function
  cf_vars <- complete_function(
    df_vars,
    n_boot = n_boot, n_perm = n_perm, alpha = alpha,
    ntree = ntree, seed = seed, n_cores = n_cores,
    verbose = verbose, always_predictors = always_predictors, which = which, plot = FALSE,
    categorical_thr = categorical_thr,
    quantitative_thr = quantitative_thr,
    importance_method = importance_method,
    prob = prob
  )

  # Generate causal graph if requested
  p1 <- NULL
  if (plot) {
    g1 <- try(plot_causal_graph_igraph(cf_vars$robust, cf_vars$binary), silent = TRUE)
    p1 <- recordPlot()
  }

  # Return results including PCA, top variables, causal scan, and plot
  list(
    pca              = pca,
    top_variables    = vars_subset,
    scan_topvars     = cf_vars,
    plot             = p1
  )
}

