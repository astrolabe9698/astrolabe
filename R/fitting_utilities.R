#' Random Forest Tuning over mtry with %IncMSE Importances
#'
#' Fits multiple Random Forests over a grid of \code{mtry} values and selects
#' the model with the lowest in-sample MSE. Returns the best model and its
#' variable importances (\%IncMSE).
#'
#' @param x A data frame or matrix of predictors (rows = samples, cols = features).
#' @param y A numeric response vector (same length as \code{nrow(x)}).
#' @param mtry_grid Integer vector of \code{mtry} values to try.
#'   If \code{NULL}, a small grid is generated automatically.
#' @param ntree Integer. Number of trees for each Random Forest (default \code{500}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{model}: the best \code{randomForest} object found;
#'   \item \code{feature_importance}: a one-column matrix with rownames = predictors
#'         and column \code{"\%IncMSE"}.
#' }
#'
#' @details
#' The score minimized is the in-sample Mean Squared Error on \code{y} vs
#' \code{predict(model, x)}. Importances are extracted from the best model via
#' \code{randomForest::importance()} and the \code{"\%IncMSE"} column is returned.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- as.data.frame(matrix(rnorm(200*5), 200, 5))
#' y <- x[[1]] * 2 + rnorm(200)
#' out <- tune_rf(x, y, mtry_grid = 1:5, ntree = 300)
#' str(out$feature_importance)
#' }
#'
#' @importFrom randomForest randomForest importance
#' @importFrom stats predict
#' @export
tune_rf <- function(x, y, mtry_grid = 1:floor(sqrt(ncol(x))), ntree = 500) {

  # If no grid for mtry is provided, create one with 5 evenly spaced values
  if (is.null(mtry_grid)) {
    mtry_grid <- unique(pmax(1, floor(seq(1, ncol(x), length.out = 5))))
  }

  # Initialize best model and best score
  best_model <- NULL
  best_score <- Inf

  # Loop over candidate mtry values
  for (m in mtry_grid) {
    # Train random forest with current mtry
    model <- randomForest::randomForest(
      x = x,
      y = y,
      ntree = ntree,
      mtry = m,
      importance = TRUE
    )

    # Predict on training data
    pred <- predict(model, x)
    # Compute mean squared error
    score <- mean((y - pred)^2)  # Mean Squared Error

    # Update best model if score improves
    if (score < best_score) {
      best_score <- score
      best_model <- model
    }
  }

  # --- Feature importance ---
  importance_matrix <- randomForest::importance(best_model)

  # Ensure importance is a proper matrix (in case of single variable)
  if (is.null(dim(importance_matrix))) {
    importance_matrix <- matrix(importance_matrix, ncol = 1)
    rownames(importance_matrix) <- names(importance_matrix)
    colnames(importance_matrix) <- "value"
  }

  # Extract importance based on %IncMSE
  importance_values <- importance_matrix[, "%IncMSE", drop = FALSE]

  # Return best model and importance values
  return(list(
    model = best_model,
    feature_importance = importance_values
  ))
}


#' Causal Entropy Scan over All Predictor→Outcome Combinations
#'
#' For each candidate outcome in \code{df}, fits a tuned Random Forest to predict
#' that outcome from the remaining variables, computes residuals, and evaluates a
#' multivariate entropy \code{entropy_nd} on \code{cbind(predictors, residual)}.
#' The final score per relation is \code{H = exp(-e_complete)}; larger is better.
#'
#' @param df Data frame containing all variables.
#' @param ntree Integer. Trees for Random Forest fitting (default \code{500}).
#' @param mtry_grid Integer vector for \code{mtry} search; default \code{1:sqrt(p)}
#'   with \code{p = ncol(df)-1}.
#' @param verbose Logical. Print per-combination diagnostics (default \code{FALSE}).
#' @param fixed_outcome Optional character scalar. If provided, only this column
#'   is used as outcome; otherwise all columns except \code{always_predictors}
#'   are considered outcomes in turn.
#' @param always_predictors Optional character vector of predictors that must
#'   always be available; excluded from the outcome set when scanning.
#' @param as_factor Optional character vector of column names to coerce to factor
#'   before fitting (also appended to \code{always_predictors}).
#' @param drill Logical. If TRUE, stops the scan after the first combination (useful for drill-down).
#'
#' @return An invisible list with:
#' \itemize{
#'   \item \code{entropy}: named numeric vector of \code{H} scores with names
#'         like \code{"X1 + X2 → Y"};
#'   \item \code{importances}: list of per-relation named lists of RF importances
#'         (one numeric value per predictor).
#' }
#'
#' @details
#' For each \code{outcome_col}, predictors are \code{setdiff(colnames(df), outcome_col)}.
#' A best RF is selected by \code{\link{tune_rf}}. Residuals are computed as
#' \code{outcome - predict(RF)}. The entropy \code{e_complete} is calculated by
#' \code{entropy_nd} on the numeric subset of \code{cbind(predictors, residual)} with
#' \code{normalize = "divide"}. The returned score is \code{H = exp(-e_complete)}.
#'
#' @examples
#' \dontrun{
#' # df must contain only the variables to scan; factors are allowed.
#' res <- causal_entropy_combinations(df, ntree = 300, verbose = TRUE)
#' head(res$entropy[order(res$entropy, decreasing = TRUE)])
#' }
#'
#' @seealso \code{\link{tune_rf}}, \code{entropy_nd}
#' @importFrom stats predict sd
#' @export
causal_entropy_combinations <- function(df,
                                        ntree = 500,
                                        mtry_grid = 1:sqrt(ncol(df)-1),
                                        verbose = FALSE,
                                        fixed_outcome = NULL,
                                        always_predictors = NULL,
                                        as_factor = NULL,
                                        drill = FALSE) {
  if (!is.null(as_factor)){
    df[,as_factor] <- as.factor(df[,as_factor])
    always_predictors <- c(always_predictors,as_factor)
  }
  H_mixed_vector <- NULL
  list_of_importances <- list()
  if (!is.null(always_predictors)) {
    always_predictors <- intersect(always_predictors, colnames(df))
  }

  outcomes <- if (is.null(fixed_outcome)) {
    setdiff(colnames(df), always_predictors)
  } else {
    fixed_outcome
  }

  to_remove <-list()
  for (outcome_col in outcomes) {
    if_res <- TRUE
    to_remove <- list()
    predictor_cols <- setdiff(colnames(df), outcome_col)
    while (TRUE){
        predictor_cols <- setdiff(predictor_cols, unlist(to_remove))
    if (length(predictor_cols) < 1) {
      if_res <- FALSE
      break;
    }
    df_predictors <- df[, predictor_cols, drop = FALSE]
    outcome_vector <- df[[outcome_col]]
    cols <- c(predictor_cols,outcome_col)

    rf <- tune_rf(df_predictors, outcome_vector, mtry_grid, ntree)
    model_rf <- rf$model
    importance_rf <- as.list(rf$feature_importance[,1])
    names(importance_rf) <- colnames(df_predictors)


    predictions <- predict(model_rf, df_predictors)
    residuals <- outcome_vector - predictions

    entropy_data <- cbind(df_predictors, residual = residuals)
    entropy_data <- entropy_data[,sapply(colnames(entropy_data),function(x){!is.factor(entropy_data[[x]])})]
    e_complete <- entropy_nd(entropy_data, normalize="divide")
    importance_rf <- unlist(importance_rf)

    H_mixed <-  exp(- e_complete)


    if (verbose) {
      cat('\n')
      cat(" Combination:", paste(predictor_cols, collapse = " + "), "→", outcome_col, "\n")
      cat(" e_complete:", e_complete, "\n")
      cat(" H_mixed   :", H_mixed, "\n")
      cat(" Mean_imp: ", mean(unlist(importance_rf)),"\n")
      cat(" Importances of the variables: ",paste(names(importance_rf),': ',round(importance_rf,2), collapse = '   ', sep = ''))
      cat('\n')
    }

    if (drill) break;
    if(all(importance_rf>0)) break;
    to_remove[[outcome_col]] <- names(importance_rf)[importance_rf<0]
      }

    if(if_res){
    model_testing <- paste(paste(predictor_cols, collapse = " + "),outcome_col,sep = " → ")
    list_of_importances[[model_testing]] <- importance_rf
    H_mixed_vector[[model_testing]] <- H_mixed
    } else {
      next
    }
  }
  invisible(list(entropy = unlist(H_mixed_vector),
                 importances = list_of_importances)
  )
}
