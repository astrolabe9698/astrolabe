#' Draw a small DAG with smart padding and optional curved edges
#'
#' @description
#' Plots a directed acyclic graph (DAG) from an edge list using **ggraph**,
#' with automatic axis limits padding and optional per-edge curvature.
#' Edges can be color-coded via a `custom_color` column (e.g., `"pairwise"`,
#' `"complex"`, `"both"`), and a subset of edges can be drawn as arcs.
#'
#' @param edges A `data.frame` with at least two columns `from` and `to`
#'   (character or factor) describing directed edges. If present, a column
#'   `custom_color` is used to map edge colors via a manual palette.
#' @param curved One of:
#'   * `NULL` (default): no curved edges;
#'   * a logical vector of length `nrow(edges)` marking which edges are curved;
#'   * a character vector of keys `"X->Y"` selecting edges to curve;
#'   * a `data.frame` with columns `from` and `to` selecting edges to curve.
#' @param layout A layout name passed to `ggraph::create_layout()` / `ggraph()`
#'   (e.g., `"auto"`, `"kk"`, `"fr"`, `"sugiyama"`, `"linear"`, ...).
#' @param pad Numeric padding added around the computed x/y ranges to
#'   prevent clipping (applied symmetrically on both axes).
#' @param arrow_len_pt Arrow length (in points) for directed edges.
#' @param end_cap_mm End cap radius (in millimeters) for edge arrows.
#' @param linewidth Edge line width.
#' @param node_size Node point size.
#' @param node_stroke Node point stroke width.
#' @param strength_curved Curvature strength for curved edges passed to
#'   `ggraph::geom_edge_arc2()` (non-curved edges use 0).
#'
#' @details
#' The function computes an automatic bounding box based on the chosen
#' `layout` and expands it by `pad` on both axes to reduce clipping and
#' keep the graph compact but readable.
#' If `edges$custom_color` exists, it is mapped with a fixed manual scale:
#' `"pairwise"` → grey, `"complex"` → black, `"both"` → greenish.
#'
#' The `curved` argument supports multiple convenient notations. When a
#' character vector is supplied, edges are identified by the key
#' `paste(from, to, sep = "->")`.
#'
#' @return A **ggplot** object.
#'
#' @examples
#' # Minimal example
#' edges <- data.frame(
#'   from = c("X","Z","Z","A"),
#'   to   = c("Y","Y","X","B"),
#'   custom_color = c("pairwise","both","complex","complex")
#' )
#'
#' # Curving a specific edge by key:
#' p1 <- draw_dag(edges, curved = "Z->X", layout = "kk")
#' # Curving by logical vector:
#' p2 <- draw_dag(edges, curved = c(FALSE, TRUE, FALSE, TRUE), layout = "fr")
#' # Curving via data.frame(from, to):
#' sel <- data.frame(from = "Z", to = "X")
#' p3 <- draw_dag(edges, curved = sel, layout = "sugiyama")
#'
#' # Print one:
#' # print(p1)
#'
#' @seealso [ggraph::ggraph()], [ggraph::geom_edge_arc2()], [igraph::graph_from_data_frame()]
#' @export
draw_dag <- function(
    edges,
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
  stopifnot(is.data.frame(edges), all(c("from","to") %in% names(edges)))

  # Mark which edges should be curved
  key <- paste(edges$from, edges$to, sep = "->")
  curved_flag <- rep(FALSE, nrow(edges))
  if (!is.null(curved)) {
    if (is.logical(curved) && length(curved) == nrow(edges)) {
      curved_flag <- curved
    } else if (is.character(curved)) {
      curved_flag <- key %in% curved
    } else if (is.data.frame(curved) && all(c("from","to") %in% names(curved))) {
      curved_key <- paste(curved$from, curved$to, sep = "->")
      curved_flag <- key %in% curved_key
    } else {
      stop("Argument 'curved' not recognized: use logical, character 'X->Y', or a data.frame with from/to.")
    }
  }
  edges$curved <- curved_flag
  strength_val <- ifelse(edges$curved, strength_curved, 0)

  # Build the graph with the edge attribute included
  g <- graph_from_data_frame(edges, directed = TRUE)

  # Tight auto limits based on the chosen layout
  auto_limits <- function(g, layout = "kk", pad = 0.4) {
    lay <- ggraph::create_layout(g, layout = layout)
    xr <- range(lay$x, na.rm = TRUE); yr <- range(lay$y, na.rm = TRUE)
    dx <- diff(xr); dy <- diff(yr)
    px <- ifelse(dx == 0, 1, dx * pad)
    py <- ifelse(dy == 0, 1, dy * pad)
    list(xlim = xr + c(-px, px), ylim = yr + c(-py, py))
  }
  lims <- auto_limits(g, layout = layout, pad = pad)

  p <- ggraph(g, layout = layout) +
    geom_edge_arc2(
      aes(edge_colour = custom_color),
      strength = strength_val,
      arrow   = grid::arrow(type = "closed", length = unit(arrow_len_pt, "pt")),
      end_cap = circle(end_cap_mm, "mm"),
      linewidth = linewidth) +
    scale_edge_color_manual(name = 'Type', values = c('pairwise' = 'grey65','complex' = 'black','both' = '#32a852'))

  p + geom_node_point(size = node_size, shape = 21, fill = "#bae3c5",
                      color = "black", stroke = node_stroke) +
    geom_node_text(aes(label = name), size = 5, color = "black",
                   vjust = 0.5, fontface = "bold") +
    theme_void() +
    coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    theme(
      legend.position   = "top",
      legend.box        = "vertical",
      legend.box.just   = "center",
      legend.box.margin = margin(t = 8, r = 0, b = 0, l = 0),  # spacing between panel and legend
      plot.margin       = margin(t = 16, r = 16, b = 12, l = 16),
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 12)
    )
}
