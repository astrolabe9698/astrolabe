#' Recursive Drill-Down over Candidate Predictor Sets
#'
#' Starting from a current best relation, iteratively explores reduced predictor
#' subsets (removing only those marked as \emph{removable}) to seek higher-entropy
#' combinations for a fixed outcome. At each step it re-fits the causal model and
#' keeps track of entropy gains and importances across layers.
#'
#' @param df Data frame with predictors and the fixed outcome column.
#' @param verdict Named numeric vector of current best relations with their entropy
#'   (e.g., \code{"X + Y → Z" = 0.42}). Passed forward and extended during drilling.
#' @param importances List accumulating variable-importance structures per layer.
#'   Typically starts from the first-layer importances returned by
#'   \code{causal_entropy_combinations()}.
#' @param layer_level Integer index of the current drill layer (e.g., 1 for first drill step).
#' @param removable_predictors Character vector of predictors deemed removable at this step.
#' @param ntree Integer. Number of trees for the underlying Random Forest (default \code{500}).
#' @param verbose Logical. Print step-by-step diagnostics (default \code{FALSE}).
#' @param fixed_outcome Character scalar. The outcome variable name kept fixed during drilling.
#' @param always_predictors Optional character vector of predictors that must always be included.
#' @param quantitative_thr Numeric. Threshold for quantitative variables (importance-based pruning).
#' @param categorical_thr Numeric. Threshold for categorical variables (importance-based pruning).
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"}.
#'   Strategy used by \code{evaluate_importance()} to decide removability.
#' @param prob Numeric. Probability cutoff for the neural network inside
#'   \code{evaluate_importance()} when \code{importance_method = "net_clust"} (default \code{0.75}).
#'
#' @return A (possibly nested) list with elements:
#' \itemize{
#'   \item \code{statistic}: numeric scalar, entropy of the winning combination at this layer;
#'   \item \code{decision}: character, relation label of the winner (e.g., \code{"X + Y → Z"});
#'   \item \code{data}: data frame of the winner's predictors + outcome;
#'   \item \code{verdict}: named numeric vector of surviving candidates and entropies;
#'   \item \code{importances}: list accumulating per-layer importances;
#'   \item \code{last_layer}: integer, index of the last processed layer;
#'   \item \code{drill_down}: recursive result for the next layer (or \code{NULL} if no further drilling).
#' }
#'
#' @details
#' For each removable predictor, the function fits \code{causal_entropy_combinations()}
#' on the remaining set and compares entropies, keeping the best. Importances for the
#' winning subset are then re-evaluated via \code{evaluate_importance()} to decide the
#' next removable set, and the procedure recurses until stopping criteria are met.
#'
#' @examples
#' \dontrun{
#' res <- drill_down_scan(
#'   df = mydata,
#'   verdict = c("A + B → Y" = 0.37),
#'   importances = list(n1th_layer = my_first_layer_imps),
#'   layer_level = 1,
#'   removable_predictors = c("A"),
#'   ntree = 500,
#'   verbose = TRUE,
#'   fixed_outcome = "Y",
#'   always_predictors = NULL,
#'   quantitative_thr = 40,
#'   categorical_thr = 35,
#'   importance_method = "neg_exp",
#'   prob = 0.75
#' )
#' }
#'
#' @seealso \code{\link{scan_all_outcomes_complete}}, \code{\link{evaluate_importance}},
#'   \code{\link{causal_entropy_combinations}}
#' @export
drill_down_scan <- function(df, verdict, importances, layer_level, removable_predictors,
                            ntree = 500,
                            verbose = FALSE,
                            fixed_outcome,
                            always_predictors=NULL,
                            quantitative_thr = 40,
                            categorical_thr = 35,
                            importance_method = c("fixed","neg_exp","net_clust"),
                            prob = 0.75) {

  # Select predictors excluding the fixed outcome
  predictors <- setdiff(colnames(df), fixed_outcome)

  # If there are fewer than 2 predictors, stop
  if (length(predictors) < 2) {
    if (verbose) cat("\033[90m\n🔚 Not enough tools to dig — only one predictor left.\n\033[0m")
    return(NULL)
  }

  # Initialize lists to store entropies, subsets, and importances
  H_list <- list()
  subsets_list <- list()
  temp_importances <- list()

  # If all predictors are removable, stop the drill
  if(base::setequal(predictors,removable_predictors)){
    if (verbose) cat("   ❌ All predictors deemed removable — stopping drill-down.\n\n")
    return(NULL)
  }

  # Iterate over each predictor to remove one at a time
  for (to_remove in predictors) {
    if (to_remove %in% removable_predictors) {
      if (verbose) cat("\n🟫🟫🟫⬇️ Digging by removing '", to_remove, "'...\n", sep = "")
      # Build reduced dataset without the removable predictor
      subset_vars <- setdiff(predictors, to_remove)
      if (length(subset_vars) < 1) next
      data_red <- df[, c(subset_vars, fixed_outcome)]
      key <- paste(paste(subset_vars, collapse = " + "), "→", fixed_outcome)
      subsets_list[[key]] <- data_red
      # Fit causal entropy model on reduced data
      fit <- causal_entropy_combinations(data_red, ntree = ntree, fixed_outcome = fixed_outcome,always_predictors = always_predictors)
      H_obs <- fit$entropy
      H_list[[key]] <- H_obs
      # Save importance of predictors for this subset
      temp_importances[[key]] <- fit$importances[[key]]
    } else {
      # Skip predictors that are "too precious" to remove
      if (verbose) cat("\033[90m\n⛓️  Holding '", to_remove, "' — too precious to remove now.\n\033[0m", sep = "")
    }
  }

  # Save importances for the new layer
  importances[[paste0('n', layer_level + 1, 'th_layer')]] <- temp_importances

  # If no entropy values were computed, stop
  if (length(H_list) == 0) {
    if (verbose) cat("\033[90m\n🕳️ No valid path found at this layer.\n\033[0m")
    return(NULL)
  }

  # Combine entropies into matrix form and find the maximum
  H_obs <- do.call(rbind, H_list)
  colnames(H_obs) <- NULL
  H_obs <- t(H_obs)

  # Select the best model (highest entropy)
  winner_name <- H_obs[which.max(H_obs)]
  names(winner_name) <- colnames(H_obs)[which.max(H_obs)]
  verdict_updated <- c(verdict, winner_name)

  # Retrieve the data corresponding to the winner model
  winner_data <- subsets_list[[names(winner_name)]]
  if (verbose) {
    cat("\033[96m\n🏆 Best combination: ", names(winner_name), " | H = ", winner_name, "\033[0m\n", sep = "")
  }

  # Store results of this layer
  res <- list(
    statistic = winner_name,
    decision = names(winner_name),
    data = winner_data,
    verdict = verdict_updated,
    importances = importances,
    last_layer = layer_level + 1
  )

  # Extract importance values for predictors of the winning model
  imp_vec <- unlist(temp_importances[[names(winner_name)]],use.names = TRUE)

  # Determine which predictors are removable based on importance evaluation
  imp_removable <- evaluate_importance(df = df,
                                       imp_vec = imp_vec,
                                       predictors = setdiff(colnames(winner_data),fixed_outcome),
                                       importance_method = importance_method,
                                       quantitative_thr = quantitative_thr,
                                       categorical_thr = categorical_thr,
                                       verbose = verbose,
                                       prob = prob
  )

  removable_predictors <- names(imp_removable[imp_removable])

  # Recursive call to drill further into the reduced dataset
  res$drill_down <- drill_down_scan(
    res$data,
    verdict = res$verdict,
    importances = res$importances,
    layer_level = res$last_layer,
    removable_predictors = removable_predictors,
    ntree = ntree,
    verbose = verbose,
    always_predictors = always_predictors,
    fixed_outcome = fixed_outcome,
    quantitative_thr = quantitative_thr,
    categorical_thr = categorical_thr,
    importance_method = importance_method,
    prob = prob
  )
  invisible(res)
}





#' Full Multi-Outcome Causal Scan with Drill-Down
#'
#' Fits \code{causal_entropy_combinations()} on all outcomes in \code{df},
#' ranks candidate relations by entropy, applies a drill-down search on each top
#' candidate (recursively removing low-importance predictors), and finally removes
#' bidirectional duplicates by keeping the higher-entropy direction.
#'
#' @param df Data frame containing all variables to scan (predictors and outcomes).
#' @param ntree Integer. Number of trees for Random Forest (default \code{500}).
#' @param verbose Logical. Print progress and diagnostics (default \code{FALSE}).
#' @param always_predictors Optional character vector of predictors that must always be included.
#' @param seed Optional integer random seed for reproducibility.
#' @param categorical_thr Numeric. Threshold for categorical predictors in
#'   \code{evaluate_importance()} (default \code{35}).
#' @param quantitative_thr Numeric. Threshold for quantitative predictors in
#'   \code{evaluate_importance()} (default \code{40}).
#' @param importance_method One of \code{"fixed"}, \code{"neg_exp"}, \code{"net_clust"}.
#'   Strategy to decide removability in \code{evaluate_importance()}.
#' @param prob Numeric. Probability cutoff for the neural network used when
#'   \code{importance_method = "net_clust"} (default \code{0.75}).
#'
#' @return A list of results (one per retained relation) where each element contains:
#' \itemize{
#'   \item \code{outer_layer_statistics}: named numeric vector of first-layer candidates and entropies;
#'   \item \code{outer_layer_decision}: character, root relation label examined;
#'   \item \code{drill_down_statistic}: named numeric scalar, final winning relation and its entropy;
#'   \item \code{drill_down_decision}: character, final winning relation label (after drilling);
#'   \item \code{verdict}: named numeric vector of surviving candidates at the final depth;
#'   \item \code{importances}: list of per-layer importances accumulated during drilling.
#' }
#' Returns \code{invisible(NULL)} if no candidates are found.
#'
#' @details
#' Pipeline:
#' \enumerate{
#'   \item Fit \code{causal_entropy_combinations(df, ...)} to obtain entropies and first-layer importances;
#'   \item Rank candidates by entropy and, for each, compute removable predictors via
#'         \code{evaluate_importance()};
#'   \item Run \code{drill_down_scan()} to explore reduced subsets until convergence;
#'   \item Remove bidirectional duplicates by keeping the higher-entropy direction.
#' }
#'
#' @examples
#' \dontrun{
#' results <- scan_all_outcomes_complete(
#'   df = mydata,
#'   ntree = 500,
#'   verbose = TRUE,
#'   always_predictors = NULL,
#'   seed = 123,
#'   categorical_thr = 35,
#'   quantitative_thr = 40,
#'   importance_method = "neg_exp",
#'   prob = 0.75
#' )
#' }
#'
#' @seealso \code{\link{drill_down_scan}}, \code{\link{evaluate_importance}},
#'   \code{\link{causal_entropy_combinations}}
#' @export
scan_all_outcomes_complete <- function(df,
                                       ntree = 500,
                                       verbose = FALSE, always_predictors = NULL, seed = NULL,
                                       categorical_thr = 35,
                                       quantitative_thr = 40,
                                       importance_method = c("fixed","neg_exp","net_clust"),
                                       prob = 0.75) {

  # Select importance method
  importance_method <- match.arg(importance_method,
                                 c("fixed","neg_exp","net_clust"))
  flag <- FALSE

  # Set seed if provided
  if(!is.null(seed)) set.seed(seed)

  # Fit entropy models for all predictor → outcome combinations
  fit <- causal_entropy_combinations(
    df,
    ntree = ntree,
    always_predictors = always_predictors
  ) #fitting to find entropies of all relationships

  H_obs <- fit$entropy #extract entropies
  if (length(H_obs) == 0) return(NULL) #quit function if there's not entropy at all

  imp_all <- list()
  imp_all[['n1th_layer']] <- fit$importances #save importances of the variables in the first layer
  ranked <- sort(H_obs, decreasing = TRUE) #rank the model basing on the entropy measure
  results <- list() #initialize list to save results

  # Verbose output: start scan
  if (verbose) {
    cat("\n────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n")
    cat("🔍 Beginning full multi-outcome scan...\n")
    cat("   • Candidates found:", length(ranked), "\n")
    cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n\n")
  }

  i <- 1
  # Loop over all candidate relationships
  for (name in names(ranked)) { #loop over the relationships found
    if (verbose) {
      cat(sprintf("➡️  [%d/%d] ROOT: %s\n", i, length(ranked), name))
      cat(sprintf("   • Initial entropy (H): %.4f\n", ranked[name]))
    }

    #extract predictors and outcomes
    parts <- strsplit(name, " → ")[[1]]
    predictors <- strsplit(parts[1], " \\+ ")[[1]]
    outcome <- parts[2]
    df_sub <- df[, c(predictors, outcome)] #take just the variables on interest (all here)

    # Get importance values of predictors for this model
    imp_vec <- unlist(imp_all[['n1th_layer']][[name]]) #extract the importances of each variables of the studied model
    names(imp_vec) <- predictors #force the name of the predictors

    # Determine removable predictors
    imp_removable <- evaluate_importance(df = df, imp_vec = imp_vec,
                                         predictors = predictors,
                                         importance_method = importance_method,
                                         quantitative_thr = quantitative_thr,
                                         categorical_thr = categorical_thr,
                                         verbose = verbose,
                                         prob = prob
    ) #evaluete importance in order to choose the less important variables

    if (verbose) {
      cat("   • Predictors:", paste(predictors, collapse = ", "), "\n")
      if (any(imp_removable)) {
        cat("   • Removable predictors:",
            paste(sprintf("%s (imp=%.3f)", names(imp_removable)[imp_removable], imp_vec[imp_removable]),
                  collapse = ", "), "\n")
      } else {
        cat("   • Removable predictors: none\n")
      }
    }

    # If all predictors are removable, skip this path
    if (sum(imp_removable) == length(predictors)) {
      if (verbose) cat("   ❌ All predictors deemed removable — skipping drill-down.\n\n")
      drill_results <- NULL
      next
    }

    drill_results <- NULL #initialize the drill
    if (length(predictors) > 1) {
      if (verbose) {
        if (any(imp_removable)) {
          cat("   ▸ Branching via 'removable' for: ",
              paste(sprintf("%s (imp=%.3f)", names(imp_removable)[imp_removable], imp_vec[imp_removable]),
                    collapse = ", "), "\n", sep = "")
        }
        cat("   ▸ Starting drill-down...\n")
      }
      #start the drill iteration
      drill_results <- drill_down_scan(
        df_sub,
        verdict = ranked[name],
        importances = imp_all,
        removable_predictors = names(imp_removable[imp_removable]),
        layer_level = 1,
        ntree = ntree,
        verbose = verbose,
        fixed_outcome = outcome,
        always_predictors = always_predictors,
        quantitative_thr = quantitative_thr,
        categorical_thr = categorical_thr,
        importance_method = importance_method,
        prob = prob
      )
    }

    # Dive into the deeper results (navigate recursive structure)
    diving_into <- drill_results
    depth <- 1
    while (!is.null(diving_into$drill_down)) {
      diving_into <- diving_into$drill_down
      depth <- depth + 1
    }

    larger <- list()
    verdict <- if (!is.null(diving_into)) diving_into$verdict else ranked[name]
    if (!is.null(diving_into)) {
      if (verbose) cat(sprintf("   ▸ Drill-down completed. Depth reached: %d\n", depth))
      # Loop over importances per layer and prune if values are too low
      for (layer in names(diving_into$importances)) {
        q=1
        for (k in names(verdict)) {
          larger[[paste(layer, q)]] <- diving_into$importances[[layer]][[k]]
          removed_candidate_name <- k
          imp_vec_layer <- unlist(larger[[paste(layer, q)]])
          under5_idx <- which(imp_vec_layer < 5)
          if (length(under5_idx) > 0 && length(verdict) > 0) {
            if (length(verdict) > 0) {
              current_best_name <- names(verdict[find_maximum(verdict)])
            } else {
              current_best_name <- NA
              break
            }
            verdict <- verdict[names(verdict) != k]

            # Verbose output if best candidate removed
            if (verbose && !is.na(current_best_name) && removed_candidate_name == current_best_name) {
              off_vars <- names(imp_vec_layer)[under5_idx]
              min_imp <- min(imp_vec_layer, na.rm = TRUE)
              cat(sprintf("   🧹 Discarded TOP candidate (post-drill, prune < 5) '%s' at layer '%s': %s (min=%.3f)\n",
                          removed_candidate_name, layer, paste(off_vars, collapse = ", "), min_imp))
              flag <- TRUE
            }
          }
          q <- q+1
        }
      }
    } else {
      # If no drill-down structure exists, prune at first layer
      if (verbose) cat("   ▸ No drill-down structure; checking first-layer importances (prune < 5).\n")
      if (any(imp_vec < 5)) {
        if (verbose) {
          off_vars <- names(which(imp_vec < 5))
          min_imp <- min(imp_vec, na.rm = TRUE)
          cat(sprintf("   🧹 Excluded root candidate (first-layer prune < 5): %s (min=%.3f)\n",
                      paste(off_vars, collapse = ", "), min_imp))
          cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n\n")
        }
        next
      }
    }

    # Extract the final decision (best model after pruning and drill-down)
    final_decision <- if (length(verdict)!=0) verdict[find_maximum(verdict)] else NULL
    final_name <- if (length(verdict)!=0) names(final_decision) else NULL
    final_H <- if (length(verdict[final_name])!=0) verdict[final_name] else NULL

    if (verbose) {
      if(!flag) cat(sprintf("🏁 Final decision: %s\n",ifelse(length(final_name) > 0, final_name, "❌ None")))
      cat(sprintf("   • Final entropy (H): %.4f\n",
                  ifelse(!is.na(final_H), final_H, NA)))
      cat(sprintf("   • Depth reached: %d layer(s)\n", depth))
      cat("────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n\n")
    }

    # Save results if decision is valid
    if(!is.null(final_decision)){
      results[[paste0("res_", i)]] <- list(
        outer_layer_statistics = ranked,
        outer_layer_decision = name,
        drill_down_statistic = final_decision,
        drill_down_decision = final_name,
        verdict = verdict,
        importances = imp_all
      )
      i <- i + 1
    }
  }

  # --- Check bidirectional duplicates (avoid keeping both directions) ---
  list_rel <- lapply(results, function(res) {
    if (!is.null(res$drill_down_decision)) {
      res$drill_down_decision
    } else {
      res$outer_layer_decision
    }
  })

  list_pred_out <- lapply(list_rel, function(rel) {
    parts <- strsplit(rel, " → ")[[1]]
    predictors <- trimws(strsplit(parts[1], "\\+")[[1]])
    out <- trimws(parts[2])
    list(p = predictors, o = out, rel = rel)
  })

  keep <- rep(TRUE, length(list_pred_out))

  if (verbose) {
    cat("🔎 Post-processing: bidirectionality check\n")
  }

  # Compare pairs of relationships for bidirectionality
  for (i in seq_along(list_pred_out)) {
    if (!keep[i]) next

    rel_i <- list_pred_out[[i]]
    preds_i <- rel_i$p
    out_i <- rel_i$o

    for (j in seq_along(list_pred_out)) {
      if (i == j || !keep[j]) next

      rel_j <- list_pred_out[[j]]
      preds_j <- rel_j$p
      out_j <- rel_j$o

      # If outcome of one is predictor of the other, check which to keep
      if ((out_i %in% preds_j) && (out_j %in% preds_i)) {
        safe_get <- function(verdict, rel) {
          if (!is.null(verdict) && rel %in% names(verdict)) {
            return(verdict[[rel]])
          } else {
            return(NULL)
          }
        }

        H_i <- safe_get(results[[i]]$verdict, list_rel[[i]])
        H_j <- safe_get(results[[j]]$verdict, list_rel[[j]])

        if (is.null(H_i) || is.null(H_j)) next

        # Keep the relation with higher entropy
        if (H_i >= H_j) {
          keep[j] <- FALSE
          if (verbose) {
            cat("   ⚖️  Bidirectional pair detected:\n")
            cat(sprintf("      • Removing: %s (H = %.4f)\n", list_rel[[j]], H_j))
            cat(sprintf("      • Keeping : %s (H = %.4f)\n", list_rel[[i]], H_i))
          }
        } else {
          keep[i] <- FALSE
          if (verbose) {
            cat("   ⚖️  Bidirectional pair detected:\n")
            cat(sprintf("      • Removing: %s (H = %.4f)\n", list_rel[[i]], H_i))
            cat(sprintf("      • Keeping : %s (H = %.4f)\n", list_rel[[j]], H_j))
          }
        }
      }
    }
  }

  # Keep only non-duplicate results
  results <- results[keep]
  if (verbose) cat("\n✅ Multi-outcome scan complete!\n")

  # Verbose: print final relationships
  if(!is.null(results)) {
    if (verbose) cat("\n FINAL CAUSAL RELATIONSHIPS:\n")
    for (i in seq_along(results)) {
      decision <- results[[i]]$drill_down_decision
      if (is.null(decision)) next

      parts <- strsplit(decision, " → ")[[1]]
      predictors <- trimws(strsplit(parts[1], "\\+")[[1]])
      outcome <- trimws(parts[2])
      all_vars <- c(predictors, outcome)
      rel_name <- paste(predictors, collapse = " + ")
      full_relation <- paste(rel_name, outcome, sep = " → ")

      if (verbose) {
        cat(sprintf("\n   🔮 %s", full_relation))
      }
    }
    if(verbose) cat("\n")
  }
  invisible(results)
}




