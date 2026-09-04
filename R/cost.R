# bioIOT cost construction, standardization and likelihood terms

#' Linear feature cost
#'
#' \eqn{C_{ij} = -\sum_k \phi_{ijk}\,\theta_k}: the cost of moving from
#' state \eqn{i} to state \eqn{j} under feature weights \eqn{\theta}.
#'
#' @param phi (K, K, F) feature array; a (K, K) matrix is promoted to a
#'   single feature.
#' @param theta Numeric length F (a scalar is recycled for 2-D \code{phi}).
#'
#' @return (K, K) cost matrix.
#' @examples
#' phi <- array(rnorm(3 * 3 * 2), c(3, 3, 2))
#' C <- make_cost(phi, c(1, -1))
#' @export
make_cost <- function(phi, theta) {
  if (is.matrix(phi)) phi <- array(as.numeric(phi), c(dim(phi), 1))
  if (is.null(dim(phi)) || length(dim(phi)) != 3L) {
    stop("`phi` must be (K, K, F)", call. = FALSE)
  }
  phi_arr <- array(as.numeric(phi), dim(phi))
  theta <- as.numeric(theta)
  F <- dim(phi_arr)[3]
  if (length(theta) != F) {
    stop(sprintf("`theta` must have length %d to match `phi` features", F),
         call. = FALSE)
  }
  if (any(!is.finite(phi_arr)) || any(!is.finite(theta))) {
    stop("`phi`/`theta` contain non-finite values", call. = FALSE)
  }
  -rowSums(sweep(phi_arr, 3, theta, `*`), dims = 2)
}

#' Row-conditional cross-entropy
#'
#' Cross-entropy between the observed row-conditional transitions \code{T}
#' and the model \eqn{Q = P/a}. Zero-mass source states are masked out;
#' the log floor keeps gradients finite as \eqn{Q \to 0}.
#'
#' @param T (K, K) observed row-conditional transition matrix.
#' @param P (K, K) model plan.
#' @param a (K,) source masses.
#' @param eps_t Log floor.
#'
#' @return Numeric loss (lower is better).
#' @examples
#' P <- matrix(0.25, 2, 2)
#' row_ce_loss(matrix(0.5, 2, 2), P, c(0.5, 0.5))
#' @export
row_ce_loss <- function(T, P, a, eps_t = 1e-12) {
  T <- as.matrix(T)
  P <- as.matrix(P)
  if (!identical(dim(T), dim(P))) {
    stop("`T` and `P` must have identical shapes", call. = FALSE)
  }
  a <- as.numeric(a)
  Q <- row_conditional(P, a, eps_t = 1e-300)
  mask <- which(a > eps_t)
  -sum(T[mask, , drop = FALSE] * log(Q[mask, , drop = FALSE] + eps_t))
}

#' Z-score standardization of state-transition features
#'
#' Z-scores each feature across all (i, j) entries. Returns the
#' standardized array plus the (mean, sd) transform so the identical
#' transform can be re-applied to validation scenarios.
#'
#' @param phi (K, K, F) feature array (a (K, K) matrix is promoted).
#'
#' @return List with \code{phi_z} ((K, K, F) standardized array) and
#'   \code{meta} (2 x F matrix: row "mean", row "sd").
#' @examples
#' phi <- array(rnorm(3 * 3 * 2), c(3, 3, 2))
#' out <- zscore_phi(phi)
#' sd(out$phi_z)
#' @export
zscore_phi <- function(phi) {
  if (is.matrix(phi)) phi <- array(as.numeric(phi), c(dim(phi), 1))
  if (is.null(dim(phi)) || length(dim(phi)) != 3L) {
    stop("`phi` must be (K, K, F)", call. = FALSE)
  }
  dims <- dim(phi)
  p <- matrix(as.numeric(phi), nrow = dims[1] * dims[2])
  if (any(!is.finite(p))) stop("`phi` contains non-finite values", call. = FALSE)
  mu <- colMeans(p)
  sd <- apply(p, 2, function(x) sqrt(mean((x - mean(x))^2)))  # population sd, paper parity
  sd[sd < 1e-12] <- 1
  pz <- (p - mu) / sd
  list(phi_z = array(pz, dims), meta = rbind(mean = mu, sd = sd))
}
