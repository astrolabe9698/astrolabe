#' Plot Causal Graph (auto-detect robust relations and/or binary matrices)
#'
#' Builds and plots a directed causal graph by auto-detecting inputs passed via
#' \code{...}: you can provide (i) robust scan results (the list produced by
#' \code{robust_scan_all_outcomes()}), and/or (ii) a binary direction matrix (or
#' the \code{$significant} matrix from \code{cor_forest_matrix_robust_perm()}).
#' If \pkg{ggraph}/\pkg{tidygraph} are available, a polished plot is produced; otherwise
#' a base-\pkg{igraph} fallback is used. Labels are drawn \emph{inside} circular nodes.
#'
#' @param ... One or more of:
#'   \itemize{
#'     \item robust scan results (a named list of per-relation lists with
#'       fields \code{predictors}, \code{outcome}, \code{significant});
#'     \item a list from \code{cor_forest_matrix_robust_perm()} or a plain numeric
#'       matrix/data frame (interpreted as a binary adjacency where 1 = edge).
#'   }
#'   You may mix multiple robust results and/or multiple matrices; they will be merged.
#' @param base_curvature Numeric in \eqn{[0,1]}. Curvature magnitude used to separate
#'   opposite-direction edges in the fallback plot (default \code{0.22}).
#' @param arrow_size,arrow_width Numeric. Arrow size and width multipliers for the
#'   fallback base-\pkg{igraph} plot (defaults \code{1.0}, \code{1.25}).
#' @param seed Integer random seed for reproducible layouts (default \code{1}).
#'
#' @return
#'   Invisibly returns:
#'   \itemize{
#'     \item a \code{tidygraph::tbl_graph} object when using the \pkg{ggraph} pipeline, or
#'     \item an \code{igraph} object when using the fallback base plot,
#'     \item \code{NULL} if no nodes/edges could be inferred.
#'   }
#'   The function draws the plot as a side effect.
#'
#' @details
#' Inputs are merged as follows:
#' \itemize{
#'   \item \strong{Robust relations:} for each significant relation, every predictor creates a
#'         directed edge \code{predictor → outcome} tagged as \code{"complex"}.
#'   \item \strong{Binary matrix:} any nonzero entry \code{[i,j]} adds \code{i → j} tagged
#'         as \code{"binary"}.
#'   \item If both sources yield the same edge, its type becomes \code{"both"} and is rendered
#'         thicker/colored accordingly.
#' }
#' Layout choice: DAGs attempt \code{layout_with_sugiyama}; otherwise \code{kk}/\code{fr}.
#' When \pkg{ggraph} is available, node diameters are sized to fit labels.
#'
#' @examples
#' \dontrun{
#' # From robust + binary:
#' TG <- plot_causal_graph_igraph(robust_results, binary_results)
#'
#' # Only binary matrix:
#' M <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
#' M["a","b"] <- 1; M["b","c"] <- 1
#' plot_causal_graph_igraph(M)
#' }
#'
#' @importFrom igraph graph_from_data_frame simplify V edge_density as_data_frame
#' @importFrom igraph layout_with_sugiyama layout_with_kk layout_with_fr vcount ecount is_dag
#' @importFrom igraph make_empty_graph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggraph scale_edge_colour_manual scale_edge_linetype_manual scale_edge_width_identity circle
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggplot2 theme_void theme element_text ggtitle margin aes
#' @importFrom grid arrow unit textGrob gpar convertWidth convertHeight grobWidth grobHeight
#' @importFrom graphics par plot title
#' @export

plot_causal_graph_igraph <- function(...,
                                     base_curvature = 0.22,
                                     arrow_size = 1.0,
                                     arrow_width = 1.25,
                                     seed = 1) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("The 'igraph' package is required.")

  aes <- ggplot2::aes

  # Utility functions
  clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))             # Restrict values between lo and hi
  sanitize_name <- function(name) gsub("[^a-zA-Z0-9_]", "_", name) # Safe variable names
  is_robust <- function(o) is.list(o) && length(o) > 0 &&
    any(vapply(o, function(x) is.list(x) && !is.null(x$predictors) && !is.null(x$outcome), logical(1)))
  is_binary <- function(o) {
    (is.list(o) && !is.null(o$significant) && is.matrix(o$significant)) ||
      is.matrix(o) || (is.data.frame(o) && is.numeric(as.matrix(o)))
  }

  # Auto-detect robust and binary results
  args <- list(...)
  robust_list <- list(); bin_mat <- NULL
  for (o in args) {
    if (is.null(o)) next
    if (is_robust(o)) {
      robust_list <- c(robust_list, list(o))
    } else if (is_binary(o)) {
      m <- if (is.list(o) && !is.null(o$significant)) o$significant else as.matrix(o)
      if (is.null(rownames(m))) rownames(m) <- seq_len(nrow(m))
      if (is.null(colnames(m))) colnames(m) <- seq_len(ncol(m))
      if (is.null(bin_mat)) {
        bin_mat <- m
      } else {
        # Merge multiple binary matrices
        rn <- union(rownames(bin_mat), rownames(m))
        cn <- union(colnames(bin_mat), colnames(m))
        M1 <- matrix(0, length(rn), length(cn), dimnames = list(rn, cn))
        M2 <- M1
        M1[rownames(bin_mat), colnames(bin_mat)] <- bin_mat
        M2[rownames(m), colnames(m)] <- m
        bin_mat <- (M1 | M2) * 1
      }
    }
  }

  # Collect all variable names for nodes
  vars_r <- if (length(robust_list))
    unique(unlist(lapply(robust_list, function(rr)
      unlist(lapply(rr, function(x) c(x$predictors, x$outcome)))))) else character(0)
  vars_b <- if (!is.null(bin_mat)) unique(c(rownames(bin_mat), colnames(bin_mat))) else character(0)
  all_vars <- unique(c(vars_r, vars_b))
  all_vars <- all_vars[!is.na(all_vars) & nzchar(all_vars)]
  if (!length(all_vars)) { cat("❌ No nodes to plot.\n"); return(invisible(NULL)) }

  # Map original names to sanitized graph variable names
  name_map  <- setNames(paste0("var_", sanitize_name(all_vars)), all_vars)
  label_map <- setNames(names(name_map), name_map)

  # Build edges from robust results
  edges_complex <- NULL
  if (length(robust_list)) {
    tmp <- list()
    for (robust_res in robust_list) for (res in robust_res) {
      if (!is.list(res) || !isTRUE(res$significant)) next
      if (is.null(res$predictors) || is.null(res$outcome)) next
      for (pred in unique(res$predictors)) {
        if (identical(pred, res$outcome)) next
        if (!pred %in% names(name_map) || !res$outcome %in% names(name_map)) next
        tmp[[length(tmp)+1]] <- data.frame(
          from = name_map[[pred]], to = name_map[[res$outcome]],
          type = "complex", stringsAsFactors = FALSE)
      }
    }
    if (length(tmp)) edges_complex <- unique(do.call(rbind, tmp))
  }

  # Build edges from binary results
  edges_binary <- NULL
  if (!is.null(bin_mat)) {
    tmp <- list()
    for (i in seq_len(nrow(bin_mat))) {
      flab <- rownames(bin_mat)[i]; if (!flab %in% names(name_map)) next
      for (j in seq_len(ncol(bin_mat))) {
        if (i == j) next
        tlab <- colnames(bin_mat)[j]; if (!tlab %in% names(name_map)) next
        if (bin_mat[i, j] == 1)
          tmp[[length(tmp)+1]] <- data.frame(
            from = name_map[[flab]], to = name_map[[tlab]],
            type = "binary", stringsAsFactors = FALSE)
      }
    }
    if (length(tmp)) edges_binary <- unique(do.call(rbind, tmp))
  }

  # Combine edges
  pieces <- Filter(Negate(is.null), list(edges_complex, edges_binary))
  if (!length(pieces)) return(invisible(igraph::make_empty_graph(directed = TRUE)))
  edges_all <- do.call(rbind, pieces)

  # Merge duplicate edges (complex + binary → both)
  key <- paste(edges_all$from, edges_all$to, sep="→")
  final_edges <- do.call(rbind, lapply(unique(key), function(k){
    idx <- which(key == k)
    row1 <- edges_all[idx[1], ]
    types <- unique(edges_all$type[idx])
    row1$type <- if (length(types) > 1) "both" else types
    row1
  }))
  final_edges <- unique(final_edges)

  # Build igraph graph and simplify
  vertices_df <- data.frame(name = unique(c(final_edges$from, final_edges$to)))
  g <- igraph::graph_from_data_frame(final_edges, directed=TRUE, vertices=vertices_df)
  g <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
  igraph::V(g)$label <- label_map[igraph::V(g)$name]

  # Layout attempt: Sugiyama → KK → Fruchterman-Reingold fallback
  set.seed(seed)
  lay <- try(igraph::layout_with_sugiyama(g)$layout, silent=TRUE)
  if (inherits(lay, "try-error")) lay <- try(igraph::layout_with_kk(g), silent=TRUE)
  if (inherits(lay, "try-error")) lay <- igraph::layout_with_fr(g)

  # Node sizing and edge width adjustments based on density
  nV <- igraph::vcount(g)
  dens <- igraph::edge_density(g, loops = FALSE)
  lab_cex <- clamp(1.30 - 0.010*nV - 0.30*dens, 0.95, 1.35)

  nE <- igraph::ecount(g)
  base_w <- clamp(1.60 - 0.30*dens, 1.10, 1.80)
  etypes <- igraph::E(g)$type
  ewidth <- rep(base_w, nE)
  ewidth[etypes == "both"] <- ewidth[etypes == "both"] + 0.6
  ewidth <- pmax(replace(ewidth, is.na(ewidth), base_w), 0.6)

  # Curvature for bidirectional edges
  edf_base <- igraph::as_data_frame(g, what="edges")
  und_key <- ifelse(edf_base$from < edf_base$to,
                    paste0(edf_base$from,"-",edf_base$to),
                    paste0(edf_base$to,"-",edf_base$from))
  curv <- numeric(nE)
  for (uk in unique(und_key[duplicated(und_key)])) {
    idx <- which(und_key == uk)
    if (length(idx) == 2) curv[idx] <- c(-base_curvature, +base_curvature)
  }

  # If ggraph is available, use it for more aesthetic plotting
  if (requireNamespace("ggraph", quietly=TRUE) &&
      requireNamespace("tidygraph", quietly=TRUE) &&
      requireNamespace("ggplot2", quietly=TRUE) &&
      requireNamespace("grid", quietly=TRUE)) {

    # Compute node diameter based on label size
    lab_vec <- igraph::V(g)$label
    size_mm_text <- 3.6 * lab_cex
    pts <- size_mm_text * (72.27/25.4)

    width_mm <- vapply(seq_along(lab_vec), function(i){
      tg <- grid::textGrob(lab_vec[i], gp = grid::gpar(fontsize = pts))
      as.numeric(grid::convertWidth(grid::grobWidth(tg), "mm", valueOnly = TRUE))
    }, numeric(1))
    height_mm <- vapply(seq_along(lab_vec), function(i){
      tg <- grid::textGrob(lab_vec[i], gp = grid::gpar(fontsize = pts))
      as.numeric(grid::convertHeight(grid::grobHeight(tg), "mm", valueOnly = TRUE))
    }, numeric(1))

    pad_mm <- 4.0
    diam_needed <- pmax(width_mm + 2*pad_mm, height_mm + 2*pad_mm)
    node_diam_mm <- clamp(max(diam_needed), 18, 42)
    node_rad_mm  <- node_diam_mm / 2

    # Prepare tidygraph and ggraph
    edf <- edf_base
    edf$type <- factor(edf$type, levels = c("complex","binary","both"))
    edf$edge_width <- ewidth
    edf$curv <- curv
    vdf <- data.frame(name = igraph::V(g)$name, label = lab_vec, stringsAsFactors = FALSE)

    TG <- tidygraph::as_tbl_graph(igraph::graph_from_data_frame(edf, directed=TRUE, vertices=vdf))
    layout_name <- if (igraph::is_dag(g)) "sugiyama" else "kk"

    # Shorten arrows to end at node boundary
    arr_mm <- clamp(2.6 + 1.0*(0.15 - dens), 2.2, 3.8)
    end_cap  <- ggraph::circle(node_rad_mm - 1.0, "mm")
    start_cap<- ggraph::circle(0.6, "mm")

    # Build plot
    p <- ggraph::ggraph(TG, layout = layout_name) +
      ggraph::geom_edge_link(
        aes(edge_colour = type, edge_linetype = type, edge_width = edge_width,
            curvature = curv),
        arrow     = grid::arrow(length = grid::unit(arr_mm, "mm"), type = "closed"),
        end_cap   = end_cap,
        start_cap = start_cap,
        lineend   = "round",
        show.legend = TRUE
      ) +
      ggraph::geom_node_point(shape = 21, size = node_diam_mm,
                              fill = "white", colour = "gray55", stroke = 0.9) +
      ggraph::geom_node_text(aes(label = label), size = size_mm_text, fontface = 2,
                             vjust = 0.5, hjust = 0.5) +
      ggraph::scale_edge_colour_manual(values = c(complex="black", binary="gray30", both="blue3"), name = "Type") +
      ggraph::scale_edge_linetype_manual(values = c(complex="solid", binary="dashed", both="solid"), name = "Type") +
      ggraph::scale_edge_width_identity(guide = "none") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Causal graph — Complex = black-solid, Binary = dashed, Both = blue-solid") +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.02, size = 12),
                     legend.position = "top",
                     legend.justification = "left",
                     plot.margin = ggplot2::margin(6, 10, 6, 10, "pt"))

    print(p)
    return(invisible(TG))
  }

  # Fallback: base igraph plot
  pt_per_mm <- 72.27/25.4
  node_diam_mm <- clamp(max(10 + 2*4, 24), 18, 42)
  node_size_pt <- node_diam_mm * pt_per_mm
  edge_cols <- ifelse(etypes=="complex", "black", ifelse(etypes=="binary","gray30","blue3"))
  edge_lty  <- ifelse(etypes=="complex", 1, ifelse(etypes=="binary", 2, 1))

  set.seed(seed)
  lay <- try(igraph::layout_with_fr(g, niter = 2000, grid = "nogrid"), silent = TRUE)
  if (inherits(lay, "try-error")) lay <- igraph::layout_with_fr(g)

  arr_size  <- clamp(0.55 - 0.15*dens, 0.32, 0.55)
  arr_width <- clamp(1.10 - 0.20*dens, 0.92, 1.10)

  op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE)
  par(mar = c(1.2, 0.6, 2.2, 0.6))
  plot(g,
       layout = lay, rescale = TRUE, asp = 0,
       vertex.label       = igraph::V(g)$label,
       vertex.label.cex   = lab_cex,
       vertex.label.font  = 2,
       vertex.label.dist  = 0,
       vertex.color       = "white",
       vertex.frame.color = "gray55",
       vertex.size        = node_size_pt,
       vertex.label.color = "black",
       edge.color         = edge_cols,
       edge.lty           = edge_lty,
       edge.width         = ewidth,
       edge.curved        = curv,
       edge.arrow.size    = min(arr_size, arrow_size),
       edge.arrow.width   = min(arr_width, arrow_width),
       curve_multiple     = TRUE)
  title("Causal graph — Complex = black-solid, Binary = dashed, Both = blue-solid",
        line = 0.6, cex.main = 1, font.main = 2)

  invisible(g)
}




