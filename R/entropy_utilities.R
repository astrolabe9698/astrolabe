#' Multivariate Differential Entropy Estimation with the Kozachenko–Leonenko Method
#'
#' This function computes the differential entropy of a multivariate dataset
#' using the k-nearest neighbor estimator (Kozachenko–Leonenko).
#'
#' @param data A data frame or a matrix of size \eqn{n \times d},
#'   where \eqn{n} is the number of observations and \eqn{d} the dimensionality.
#' @param k Number of neighbors to consider in the estimator.
#'   Must be a positive integer, typically \eqn{k \geq 2}.
#'   Default: 10.
#' @param normalize Normalization mode of the returned entropy.
#'   Can be:
#'   \itemize{
#'     \item `"none"`: raw entropy in bits;
#'     \item `"divide"`: per-dimension entropy (bits per dimension);
#'     \item `"sqrt"`: square-root scaling of entropy (non-standard heuristic).
#'   }
#'   Default: `"none"`.
#'
#' @details
#' The Kozachenko–Leonenko estimator is based on the distance to the
#' \eqn{k}-th nearest neighbor for each point, and the volume of the unit ball in \eqn{R^d}.
#' The entropy is computed in nats and then converted to bits.
#'
#' @return A single numeric value representing the estimated entropy
#'   (possibly normalized) in bits.
#'
#' @references
#' Kozachenko, L. F., & Leonenko, N. N. (1987).
#' Sample estimate of the entropy of a random vector.
#' \emph{Problemy Peredachi Informatsii}, 23(2), 9–16.
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rnorm(100 * 3), ncol = 3)
#' entropy_nd(data, k = 5, normalize = "none")
#' entropy_nd(data, k = 5, normalize = "divide")
#'
#' @name entropy_nd
#' @export
entropy_nd <- function(data, k = 10, normalize = c("none", "divide", "sqrt")) {
  normalize <- match.arg(normalize)
  data_mat <- as.matrix(data)
  n <- nrow(data_mat)
  d <- ncol(data_mat)

  # calcolo k-NN
  knn <- FNN::get.knn(data_mat, k = k)

  # per ogni punto prendo la distanza al k-esimo vicino
  eps <- knn$nn.dist[, k]
  sum_log <- sum(log(eps + 1e-10))  # stabilità numerica

  # volume della palla unitaria in R^d
  volume_ball <- pi^(d / 2) / gamma(d / 2 + 1)

  # Kozachenko–Leonenko estimator in nats
  H_nat <- -digamma(k) + digamma(n) + log(volume_ball) + d * (sum_log / n)

  # conversione da nats a bits
  H_bits <- H_nat / log(2)

  # normalizzazione
  out <- switch(
    normalize,
    "divide" = H_bits / d,
    "sqrt"   = sqrt(H_bits),
    "none"   = H_bits
  )

  if (is.null(out)) {
    out <- H_bits^(1 / d)
  }

  return(out)
}
