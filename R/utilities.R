#' ReLU (Rectified Linear Unit) activation
#'
#' Sets all negative values to zero.
#'
#' @param x Numeric vector or matrix input.
#' @return A matrix with negative values replaced by zero.
#' @export
relu <- function(x) {
  x <- as.matrix(x)
  x[x < 0] <- 0
  x
}

#' Sigmoid activation
#'
#' Computes the sigmoid transformation.
#'
#' @param x Numeric vector or matrix input.
#' @return Transformed values in (0, 1).
#' @export
sigmoid  <- function(x) 1 / (1 + exp(-x))

#' Map name/function to activation
#'
#' Internal helper to map string names to activation functions.
#'
#' @param a Either the name of the activation function
#'   (`"relu"`, `"sigmoid"`, `"tanh"`)
#'   or a custom function.
#'
#' @return A function implementing the activation.
#' @keywords internal
.as_act <- function(a) {
  if (is.function(a)) return(a)
  switch(a,
         "relu"    = relu,
         "sigmoid" = sigmoid,
         "tanh"    = tanh,
         stop(sprintf("Unknown activation: %s", a))
  )
}

#' Manual prediction with a dense neural network
#'
#' Performs a forward pass of a feed-forward neural network
#' given a list of weights and biases, applying the specified
#' activation functions at each layer.
#'
#' @param X Numeric input matrix, with rows = observations and columns = features.
#' @param weights_dense A list of length `2*L` containing, for each layer,
#'   the weight matrix (`W`) and the bias vector (`b`), in the order `[W1, b1, W2, b2, ...]`.
#' @param activations Either:
#'   - a single activation name/function (recycled for all layers), or
#'   - a vector/list of length `L` with one activation per layer.
#'   Supported names: `"relu"`, `"sigmoid"`, `"tanh"`, or custom functions.
#' @param last_activation Optional. If provided, overrides the activation
#'   of the last layer.
#' @param verbose Logical; if `TRUE`, prints layer shapes and debugging info.
#'
#' @return A numeric matrix (or vector if 1D) with the output of the network.
#' @export
predict_manual <- function(X, weights_dense, activations = "relu", last_activation = NULL, verbose = FALSE) {
  # forza X matrice numerica
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  # numero di layer Dense
  L <- length(weights_dense) / 2
  stopifnot(L == floor(L), L >= 1)

  # normalizza le attivazioni -> lista di funzioni di lunghezza L
  if (length(activations) == 1L) {
    acts <- rep(list(activations), L)
  } else if (length(activations) == L) {
    acts <- as.list(activations)
  } else {
    stop("Length 'activations' not valid.")
  }
  if (!is.null(last_activation)) acts[[L]] <- last_activation
  acts <- lapply(acts, .as_act)

  h <- X
  for (l in seq_len(L)) {
    W <- as.matrix(weights_dense[[2*l - 1]])
    b <- as.numeric(weights_dense[[2*l]])
    storage.mode(W) <- "double"

    # mantieni h come matrice double
    h <- as.matrix(h)
    storage.mode(h) <- "double"
    stopifnot(is.matrix(h), is.matrix(W))

    if (verbose) {
      cat(sprintf("Layer %d: h %dx%d, W %dx%d, b %d\n",
                  l, nrow(h), ncol(h), nrow(W), ncol(W), length(b)))
    }

    # orientazione robusta (usa isTRUE per evitare NA)
    if (isTRUE(ncol(h) == nrow(W))) {
      z <- h %*% W
    } else if (isTRUE(ncol(h) == ncol(W))) {
      z <- h %*% t(W)
    } else {
      stop(sprintf(
        "Layer %d: mismatch dimensions -> ncol(h)=%s, nrow(W)=%s, ncol(W)=%s",
        l, paste0(ncol(h), collapse=","), paste0(nrow(W), collapse=","),
        paste0(ncol(W), collapse=",")
      ))
    }

    # somma bias per colonna
    if (!isTRUE(length(b) == ncol(z))) {
      stop(sprintf("Layer %d: bias len %d vs units %d", l, length(b), ncol(z)))
    }
    z <- sweep(z, 2, b, "+")
    z <- as.matrix(z)

    # attivazione
    h <- acts[[l]](z)
  }

  drop(h)
}

#' Find the Index of the Maximum Value
#'
#' Returns the index of the maximum element in a numeric vector.
#'
#' @param v Numeric vector.
#' @return Integer scalar: the index of the maximum value.
#'   If \code{length(v) == 1}, returns \code{1}.
#'   If \code{length(v) == 0}, raises an error.
#'
#' @examples
#' find_maximum(c(3, 5, 2))   # 2
#' find_maximum(7)            # 1
#' @export
find_maximum <- function(v) {
  if (length(v) ==1 ) {
    return(1)
  } else if (length(v) < 1 ) {
    stop("Empty vector")
  }

  return(which.max(v))
}


#' Remove Outliers from a Data Frame
#'
#' Removes all rows containing outliers in any numeric column,
#' using the interquartile range (IQR) rule.
#'
#' @param df A data frame with numeric and/or non-numeric columns.
#' @return A data frame with the same columns as \code{df}, with rows
#'   removed if any numeric column contains an outlier in that row.
#'
#' @examples
#' df <- data.frame(a = c(1,2,3,100), b = c(5,6,7,8))
#' remove_outliers(df)
#'
#' @importFrom stats quantile
#' @export
remove_outliers <- function(df) {

  # Identify numeric columns (only these are tested for outliers)
  num_cols <- sapply(df, is.numeric)

  # Column-wise outlier detector using IQR fences
  is_outlier <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1 * IQR
    upper_bound <- Q3 + 1 * IQR
    x < lower_bound | x > upper_bound
  }

  # Logical matrix/vector of outlier flags per numeric column
  outliers_logical <- sapply(df[, num_cols, drop=FALSE], is_outlier)

  # Ensure a matrix shape even if there is only one numeric column
  if (is.vector(outliers_logical)) {
    outliers_logical <- matrix(outliers_logical, ncol = 1)
  }

  # Mark rows to remove if any numeric column is outlier in that row
  rows_to_remove <- apply(outliers_logical, 1, any)

  # Keep only rows without any flagged outliers
  df_clean <- df[!rows_to_remove, ]

  return(df_clean)
}


#' Negative-Exponential Threshold Function
#'
#' Computes a decaying threshold as a function of the number of predictors (\eqn{p}),
#' optionally lower-bounded by \code{L}.
#'
#' @param p Integer. Number of predictors/dimension.
#' @param thr Numeric. Initial threshold value (at \eqn{p=1}).
#' @param k Numeric. Decay rate.
#' @param L Optional numeric. Lower bound (default \code{NULL}).
#' @return Numeric vector of thresholds.
#'
#' @examples
#' thr_exp(10, thr = 40, k = 0.05)
#' thr_exp(1:5, thr = 40, k = 0.2, L = 10)
#' @export
thr_exp <- function(p, thr, k, L = NULL) {
  val <- thr * exp(-k * (p-1))
  if (!is.null(L)) pmax(L, val) else val
}

#' Evaluate Variable Importance and Decide Removability
#'
#' Classifies predictors as removable or to be kept based on importance values,
#' using one of:
#' \itemize{
#'   \item \code{"fixed"} — fixed thresholds for quantitative vs categorical;
#'   \item \code{"net_clust"} — 1D DBSCAN clustering + neural network on high cluster(s);
#'   \item \code{"neg_exp"} — adaptive exponential-decay thresholds depending on \eqn{p}.
#' }
#'
#' @param df Data frame containing the predictors.
#' @param imp_vec Named numeric vector of importance values.
#' @param predictors Character vector of predictor names (order for the output).
#' @param importance_method One of \code{"fixed"}, \code{"net_clust"}, \code{"neg_exp"}.
#' @param quantitative_thr Numeric. Threshold for quantitative predictors (default 40).
#' @param categorical_thr Numeric. Threshold for categorical predictors (default 35).
#' @param verbose Logical. Print verbose diagnostics (default \code{FALSE}).
#' @param prob Numeric. Probability cutoff for the neural network (default 0.75).
#'
#' @return Named logical vector (aligned with \code{predictors}),
#'   where \code{TRUE} means "removable".
#'
#' @examples
#' imp <- c(x1 = 10, x2 = 50)
#' df  <- data.frame(x1 = 1:10, x2 = letters[1:10])
#' evaluate_importance(df, imp, predictors = names(imp), importance_method = "fixed")
#'
#'
#' @importFrom stats sd mad median
#' @export
evaluate_importance <- function(df,imp_vec,
                                predictors,
                                importance_method = c("fixed","net_clust","neg_exp"),
                                quantitative_thr = 40,
                                categorical_thr = 35,
                                verbose = FALSE,
                                prob = 0.75
){
  weights_dense <- get_weights_dense()
  # Pretty-printers for verbose diagnostics
  show_all <- function(tag, v){
    if (!verbose) return()
    ord <- order(v, decreasing = TRUE)
    cat(sprintf("   • %s [%d]: %s\n", tag, length(v),
                paste0(predictors[ord], " (", sprintf("%.3f", v[ord]), ")", collapse=", ")))
  }

  show_group <- function(tag, mask, v){
    if (!verbose) return()
    if (any(mask, na.rm = TRUE)) {
      cat(sprintf("   • %s [%d]: %s\n", tag, sum(mask, na.rm=TRUE),
                  paste0(predictors[mask], " (", sprintf("%.3f", v[mask]), ")", collapse=", ")))
    } else {
      cat(sprintf("   • %s: (none)\n", tag))
    }
  }

  # Normalize method selection and optionally print initial importances
  importance_method <- match.arg(importance_method, c("fixed","net_clust","neg_exp"))

  if (verbose) show_all("IMPORTANCES", imp_vec)

  # Identify categorical variables among provided predictor names
  fac <- sapply(names(imp_vec),function(x) {is.factor(df[,x])})

  # --- Method: NEG_EXP (adaptive thresholds via exponential decay) ---
  if (importance_method == "neg_exp") {
    if (any(fac)) {
      imp_removable <- rep(FALSE,length(imp_vec))
      names(imp_removable) <- names(imp_vec)
      categorical_thr_adapted <- thr_exp(length(predictors), categorical_thr, 0.025)
      quantitative_thr_adapted <- thr_exp(length(predictors), quantitative_thr, 0.025)
      imp_removable[fac] <- imp_vec[fac] < categorical_thr_adapted
      imp_removable[!fac] <- imp_vec[!fac] < quantitative_thr_adapted

      if (verbose) {
        cat(sprintf("   • New thresholds: %.4f (%s), %.4f (%s)\n",
                    categorical_thr_adapted, "categorical",quantitative_thr_adapted, "quantitative"))
      }
    } else {
      quantitative_thr_adapted <- thr_exp(length(predictors), quantitative_thr, 0.025)
      imp_removable <- imp_vec < quantitative_thr_adapted
      if (verbose) {
        cat(sprintf("   • New threshold: %.4f (%s)\n",
                    quantitative_thr_adapted, "quantitative"))
      }
    }

    # --- Method: NET_CLUST (DBSCAN clustering + neural model) ---
  } else if (importance_method == "net_clust") {

    # Edge case: only one variable → directly use neural model
    if (length(imp_vec) == 1) {
      if (verbose) cat("   • Only one variable → NEURAL NETWORK\n")
      data_tree <- data.frame(IncMSE = imp_vec, p = rep(length(predictors), length(imp_vec)), var_type = as.integer(!fac))

      # Predict keep/remove using the loaded neural network; threshold at 'prob'
      pred <- predict_manual(X  = data_tree, weights_dense = weights_dense, activations = "relu", last_activation = "sigmoid", verbose = FALSE)
      pred <- as.matrix(pred, ncol=1)
      pred <- pred > prob
      imp_removable <- !(as.vector(pred))

    } else {
      # Cluster importances in 1D using an augmented DBSCAN strategy
      cl <- dbscan_1d_augmented(
        imp_vec,
        reps_per_point = 50,
        jitter_scale = 1 / mean(diff(sort(imp_vec))),
        eps_factor = 2.5*sd(imp_vec),
        plot_result = FALSE
      )$clusters

      # Single non-noise cluster: use neural model for all variables
      if (length(unique(cl)) == 1 && unique(cl) != 0) {
        if (verbose) cat("   • DBSCAN: 1 CLUSTER only → NEURAL NETWORK on ALL\n")
        data_tree <- data.frame(IncMSE = imp_vec, p = rep(length(predictors), length(imp_vec)),var_type = as.integer(!fac))

        pred <- predict_manual(X  = data_tree, weights_dense = weights_dense, activations = "relu", last_activation = "sigmoid", verbose = FALSE)
        pred <- as.matrix(pred, ncol=1)
        pred <- pred > prob
        imp_removable <- !(as.vector(pred))

        # Exactly two clusters: remove the low-mean cluster, NN on the high-mean cluster
      } else if (length(unique(cl)) == 2) {
        if (verbose) cat("   • DBSCAN: 2 CLUSTERS detected\n")

        cl_mat <- as.data.frame(cbind(imp_vec, cl))
        colnames(cl_mat) <- c("imp_vec", "clusters")

        m <- c()
        for (ci in unique(cl_mat$clusters)) {
          m <- c(m, mean(cl_mat[cl_mat$clusters == ci, "imp_vec"]))
        }
        names(m) <- unique(cl_mat$clusters)

        if (verbose) {
          cat(sprintf("   • Cluster means: %s\n",
                      paste(names(m), sprintf("%.3f", m), sep="=", collapse=", ")))
        }

        m_low <- names(m)[which.min(m)]
        imp_removable <- imp_vec %in% cl_mat[cl_mat$clusters == m_low, "imp_vec"]
        names(imp_removable) <- names(imp_vec)

        m_high <- names(m)[which.max(m)]
        vec_high <- cl_mat[cl_mat$clusters == m_high, "imp_vec"]
        vec_high <- imp_vec[imp_vec %in% vec_high]

        if (verbose) {
          cat(sprintf("   • Cluster LOW → removed (%d vars), Cluster HIGH → NEURAL NETWORK (%d vars)\n",
                      sum(imp_removable), length(vec_high)))
        }

        data_tree <- data.frame(IncMSE = vec_high, p = rep(nrow(cl_mat), length(vec_high)),var_type = as.integer(!fac[names(vec_high)] ))

        pred <- predict_manual(X  = data_tree, weights_dense = weights_dense, activations = "relu", last_activation = "sigmoid", verbose = FALSE)
        pred <- as.matrix(pred, ncol=1)
        pred <- pred > prob
        imp_removable[names(vec_high)] <- !(as.vector(pred))

        # More than two clusters: keep only the highest-mean cluster for NN, others removable
      } else if (length(unique(cl)) > 2) {
        if (verbose) cat(sprintf("   • DBSCAN: %d CLUSTERS detected\n", length(unique(cl))))

        cl_mat <- as.data.frame(cbind(imp_vec, cl))
        colnames(cl_mat) <- c("imp_vec", "clusters")

        m <- c()
        for (ci in unique(cl_mat$clusters)) {
          m <- c(m, mean(cl_mat[cl_mat$clusters == ci, "imp_vec"]))
        }
        names(m) <- unique(cl_mat$clusters)

        if (verbose) {
          cat(sprintf("   • Cluster means: %s\n",
                      paste(names(m), sprintf("%.3f", m), sep="=", collapse=", ")))
        }

        imp_removable <- rep(TRUE, length(imp_vec))
        names(imp_removable) <- names(imp_vec)

        m_high <- names(m)[which.max(m)]
        vec_high <- cl_mat[cl_mat$clusters == m_high, "imp_vec"]
        vec_high <- imp_vec[imp_vec %in% vec_high]

        if (verbose) {
          cat(sprintf("   • Highest-mean cluster kept for NEURAL NETWORK (%d vars)\n", length(vec_high)))
        }

        data_tree <- data.frame(IncMSE = vec_high, p = rep(nrow(cl_mat), length(vec_high)),var_type = as.integer(!fac[names(vec_high)] ))

        pred <- predict_manual(X  = data_tree, weights_dense = weights_dense, activations = "relu", last_activation = "sigmoid", verbose = FALSE)
        pred <- as.matrix(pred, ncol=1)
        pred <- pred > prob
        imp_removable[names(vec_high)] <- !(as.vector(pred))

        # All points considered noise: fallback to neural model for all
      } else {
        if (verbose) cat("   • DBSCAN: all NOISE → NEURAL NETWORK on ALL\n")
        data_tree <- data.frame(IncMSE = imp_vec, p = rep(length(predictors), length(imp_vec)),var_type = as.integer(!fac))
        pred <- predict_manual(X  = data_tree, weights_dense = weights_dense, activations = "relu", last_activation = "sigmoid", verbose = FALSE)
        pred <- as.matrix(pred, ncol=1)
        pred <- pred > prob
        imp_removable <- !(as.vector(pred))
      }
    }

    # --- Method: FIXED (static thresholds per variable type) ---
  } else {
    if (any(fac)) {
      imp_removable <- rep(FALSE,length(imp_vec))
      names(imp_removable) <- names(imp_vec)

      imp_removable[fac] <- imp_vec[fac] < categorical_thr
      imp_removable[!fac] <- imp_vec[!fac] < quantitative_thr
    } else {
      imp_removable <- imp_vec < quantitative_thr
    }
  }

  # Optional verbose summary of final decision
  if (verbose) {
    show_group("REMOVED (final)", imp_removable, imp_vec)
    show_group("KEPT (final)", !imp_removable, imp_vec)
  }

  # Ensure names align with the provided predictor order
  names(imp_removable) <- predictors
  return(imp_removable)
}

#' One-Dimensional DBSCAN with Data Augmentation
#'
#' Runs DBSCAN on one-dimensional data, augmented by jittered replicates
#' around each observation to improve clustering stability.
#'
#' @param x Numeric vector of data points (no \code{NA}).
#' @param reps_per_point Integer. Number of augmented replicates per point.
#' @param jitter_scale Numeric. Scale of jitter relative to MAD (default \code{0.25}).
#' @param eps_factor Numeric. \code{eps = eps_factor * sd_jitter} (default \code{2.5}).
#' @param minPts_base Integer. Base \code{minPts} for DBSCAN (default \code{2}).
#' @param scale_data Logical. Standardize data before clustering.
#' @param seed Integer random seed (default \code{123}).
#' @param plot_result Logical. Plot the clustering result (default \code{TRUE}).
#' @param show_aug_points Logical. Show augmented points in the plot.
#' @param main Character. Plot title.
#'
#' @return List with components:
#' \item{clusters}{Integer vector of cluster assignments per original point (0 = noise).}
#' \item{noise}{Indices of noise points.}
#' \item{n_clusters}{Number of clusters detected.}
#' \item{parameters}{List of parameters used.}
#'
#' @examples
#' res <- dbscan_1d_augmented(c(1,2,3,10), plot_result = FALSE)
#' res$clusters
#'
#'
#' @importFrom stats mad median sd rnorm
#' @importFrom graphics plot points text legend abline
#' @importFrom grDevices rgb
#' @importFrom dbscan dbscan
#' @export
dbscan_1d_augmented <- function(
    x,
    reps_per_point = 60,
    jitter_scale   = 0.25,    # sd_jitter = jitter_scale * MAD(x)
    eps_factor     = 2.5,     # eps = eps_factor * sd_jitter
    minPts_base    = 2,
    scale_data     = FALSE,
    seed           = 123,
    plot_result    = TRUE,
    show_aug_points= FALSE,
    main           = "Clustering 1D con DBSCAN (augmented)"
) {

  # Basic input validation
  if (!is.numeric(x)) stop("Vector should be numeric.")
  if (anyNA(x)) stop("The vector has NAs.")
  n <- length(x); if (n < 2) stop("There must be at least two points.")

  set.seed(seed)
  x_use <- if (scale_data) as.numeric(scale(x)) else as.numeric(x)

  # Robust spread estimate (MAD); fallback to sd, then tiny epsilon if needed
  mad_x <- stats::mad(x_use, constant = 1, center = median(x_use))
  if (!is.finite(mad_x) || mad_x == 0) {
    mad_x <- stats::sd(x_use)
    if (!is.finite(mad_x) || mad_x == 0) mad_x <- 1e-6
  }
  sd_jitter <- jitter_scale * mad_x

  # --- Data augmentation: replicate each x[i] with Gaussian jitter ---
  total <- n * reps_per_point
  X_aug <- numeric(total)
  owners <- integer(total)
  idx <- 1
  for (i in seq_len(n)) {
    block <- x_use[i] + rnorm(reps_per_point, 0, sd_jitter)
    X_aug[idx:(idx + reps_per_point - 1)] <- block
    owners[idx:(idx + reps_per_point - 1)] <- i
    idx <- idx + reps_per_point
  }

  # DBSCAN parameters derived from jitter scale
  eps    <- eps_factor * sd_jitter
  minPts <- max(minPts_base, ceiling(0.6 * reps_per_point))

  # Run DBSCAN on the augmented 1D data (as a column matrix)
  db <- dbscan::dbscan(matrix(X_aug, ncol = 1), eps = eps, minPts = minPts)

  # Map augmented cluster labels back to the original n points by majority vote
  lab_orig <- integer(n)
  for (i in seq_len(n)) {
    labs_i <- db$cluster[owners == i]
    tbl <- table(labs_i[labs_i > 0])
    lab_orig[i] <- if (length(tbl)) as.integer(names(tbl)[which.max(tbl)]) else 0L
  }
  # Relabel clusters to consecutive integers (1..K), keep 0 for noise
  uniq_pos <- sort(unique(lab_orig[lab_orig > 0]))
  relabel  <- setNames(seq_along(uniq_pos), uniq_pos)
  clusters <- ifelse(lab_orig == 0, 0L, relabel[as.character(lab_orig)])

  # Pack results and the parameters actually used
  result <- list(
    clusters   = clusters,
    noise      = which(clusters == 0),
    n_clusters = length(unique(clusters[clusters > 0])),
    parameters = list(
      reps_per_point = reps_per_point,
      jitter_scale   = jitter_scale,
      sd_jitter      = sd_jitter,
      eps_factor     = eps_factor,
      eps            = eps,
      minPts         = minPts,
      scale_data     = scale_data,
      seed           = seed
    )
  )

  # Optional visualization: original points colored by cluster/noise label
  if (plot_result) {
    y0 <- rep(0, n)
    cols <- ifelse(clusters == 0, "grey50", clusters + 1L)

    plot(x, y0, type = "n",
         xlab = "Values", ylab = "", yaxt = "n", main = main)
    abline(h = 0, col = "black")

    if (show_aug_points) {
      points(
        if (scale_data) (X_aug * sd(x, na.rm=TRUE) + mean(x, na.rm=TRUE)) else X_aug,
        jitter(rep(0, length(X_aug)), amount = 0.01),
        pch = 16, cex = 0.4, col = rgb(0,0,0,0.15)
      )
    }

    points(x, y0, pch = 19, cex = 1.6, col = cols)
    text(x, y0, labels = x, pos = 3, cex = 0.8, col = cols)
    if (result$n_clusters > 0) {
      legend_labels <- c("noise (0)", paste("cluster", 1:result$n_clusters))
      legend_colors <- c("grey50", 2:(result$n_clusters+1))
    } else {
      legend_labels <- "noise (0)"
      legend_colors <- "grey50"
    }

    legend("topleft",
           legend = legend_labels,
           pch = 19,
           pt.cex = 1.0,
           cex = 0.9,
           col = legend_colors,
           bty = "n",
           y.intersp = 0.3,
           x.intersp = 0.5
    )

  }

  return(result)
}

#' Heuristic Check for Numeric Coercibility
#'
#' Tests whether a vector can be safely coerced to numeric without excessive missing values.
#'
#' @param v A vector of any type.
#' @param na_tol Numeric in \eqn{[0,1]}. Maximum tolerated fraction of \code{NA}s after coercion (default \code{0.01}).
#'
#' @return Logical. \code{TRUE} if coercion is safe, \code{FALSE} otherwise.
#'
#' @examples
#' can_coerce_numeric(c("1", "2", "3"))   # TRUE
#' can_coerce_numeric(c("1", "a", "3"))   # FALSE (di solito)
#' @export
can_coerce_numeric <- function(v, na_tol = 0.01) {
  if (is.numeric(v)) return(TRUE)
  vv <- as.character(v)
  suppressWarnings(num <- as.numeric(vv))
  na_rate <- mean(is.na(num))
  is.finite(sum(num, na.rm = TRUE)) && na_rate <= na_tol
}
