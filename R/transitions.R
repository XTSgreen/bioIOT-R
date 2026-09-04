# bioIOT transition matrices, pseudotime and state features

#' State-to-state transition matrix from a fit
#'
#' Re-solves the transport plan at the fitted \eqn{\theta} and returns the
#' row-conditional transition matrix \eqn{Q = P/a}. Works on the stored
#' scenarios of a \code{\link{fit_iot}} result (no extra arguments needed),
#' or on user-supplied \code{phi}/\code{a}/\code{b}.
#'
#' @param fit A \code{bioIOT_fit} object.
#' @param which Scenario index (default 1).
#' @param phi,a,b Optional user-supplied scenario overriding the stored one
#'   (\code{phi} must already be standardized the same way as in the fit,
#'   e.g. via \code{\link{zscore_phi}}).
#'
#' @return (K, K) row-stochastic transition matrix.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_iot_states(K = 5, seed = 1)
#' fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
#' Q <- transition_matrix(fit)
#' rowSums(Q)
#' @export
transition_matrix <- function(fit, which = 1, phi = NULL, a = NULL, b = NULL) {
  if (!inherits(fit, "bioIOT_fit")) {
    stop("`fit` must be a bioIOT_fit object (see fit_iot)", call. = FALSE)
  }
  if (is.null(phi)) {
    sc <- fit$scenarios[[which]]
    phi <- sc$phi
    a <- sc$a
    b <- sc$b
  } else {
    phi <- if (is.matrix(phi)) array(as.numeric(phi), c(dim(phi), 1)) else as.array(phi)
    a <- as.numeric(a)
    b <- as.numeric(b)
    K <- dim(phi)[1]
    a <- .prob_vector(a, K, "a")
    b <- .prob_vector(b, K, "b")
  }
  C <- make_cost(phi, fit$theta)
  P <- soft_sinkhorn(C, a, b, mu = fit$hyper$mu, eps = fit$hyper$eps,
                     iters = fit$hyper$iters, damp = fit$hyper$damp)
  Q <- row_conditional(P, a)
  rownames(Q) <- colnames(Q) <- if (!is.null(dimnames(phi)) &&
                                    !is.null(dimnames(phi)[[2]])) {
    dimnames(phi)[[2]]
  } else {
    sprintf("S%d", seq_len(nrow(Q)))
  }
  Q
}

#' Random-walk pseudotime from a transition matrix
#'
#' Expected number of random-walk steps to reach the root state under the
#' transition matrix \eqn{Q} (root absorbing). A simple, dependency-free
#' "last mile" that turns an IOT transition matrix into a per-state
#' pseudotime.
#'
#' @param Q (K, K) row-stochastic transition matrix (e.g. from
#'   \code{\link{transition_matrix}}).
#' @param root Root state: integer index or a name matching rownames of Q.
#'
#' @return Named numeric vector of expected hitting times (root = 0,
#'   increasing away from the root).
#'
#' @examples
#' Q <- matrix(c(0.5, 0.5, 0, 0, 0.2, 0.8, 0, 0, 0.1, 0.1, 0.4, 0.4, 1, 0, 0, 0), 4, 4, byrow = TRUE)
#' pseudotime_from_transition(Q, root = 1)
#' @export
pseudotime_from_transition <- function(Q, root) {
  Q <- as.matrix(Q)
  if (any(!is.finite(Q)) || any(Q < 0)) {
    stop("`Q` must be finite and non-negative", call. = FALSE)
  }
  K <- nrow(Q)
  nms <- if (!is.null(rownames(Q))) rownames(Q) else sprintf("S%d", seq_len(K))
  if (is.character(root)) {
    if (!root %in% nms) {
      stop("`root` name not found in rownames(Q)", call. = FALSE)
    }
    root <- match(root, nms)
  }
  if (!is.numeric(root) || length(root) != 1L || root < 1 || root > K ||
      root != as.integer(root)) {
    stop("`root` must be an integer index in [1, K] or a state name",
         call. = FALSE)
  }
  free <- setdiff(seq_len(K), root)
  nms <- if (!is.null(rownames(Q))) rownames(Q) else sprintf("S%d", seq_len(K))
  if (length(free) == 0L) {
    return(setNames(numeric(1), nms[root]))
  }
  # absorbing-closure: states with no outgoing mass (e.g. zero-mass source
  # rows of an IOT plan) would otherwise be transient with hitting time 1;
  # bounce them back along their inflow distribution so the walk stays
  # ergodic and their hitting time reflects being a terminal state
  zero <- rowSums(Q) < 1e-12
  if (any(zero)) {
    cs <- colSums(Q)
    for (j in which(zero)) {
      r <- Q[, j] / pmax(cs[j], 1e-12)   # inflow distribution into state j
      if (sum(r) <= 0 || !is.finite(sum(r))) r <- rep(1, K)
      Q[j, ] <- r / sum(r)
    }
  }
  Qf <- Q[free, free, drop = FALSE]
  # absorbing-chain equation (I - Qf) h = 1; Qf is raw (sub-stochastic),
  # mass flowing to the root is what makes the system well-posed
  h <- as.vector(solve(diag(length(free)) - Qf, rep(1, length(free))))
  pt <- numeric(K)
  pt[free] <- h
  pt[root] <- 0
  setNames(pt, nms)
}

#' Build IOT state-transition features from state embeddings
#'
#' Builds the (K, K, F) feature array used by the paper: per-dimension
#' "pure-column" target-state features plus a state-similarity interaction
#' block \eqn{\phi_{ij,F+1} = u_i \cdot u_j}.
#'
#' @param u (K, D) matrix of state feature vectors (e.g. cluster centroids
#'   in a reduced space). Column names become feature labels.
#'
#' @return (K, K, D + 1) feature array with dimnames
#'   \code{list(NULL, NULL, c(colnames(u), "sim"))}.
#'
#' @examples
#' u <- matrix(rnorm(6 * 2), 6, 2, dimnames = list(NULL, c("PC1", "PC2")))
#' phi <- build_state_features(u)
#' dim(phi)
#' @export
build_state_features <- function(u) {
  u <- as.matrix(u)
  if (any(!is.finite(u))) stop("`u` contains non-finite values", call. = FALSE)
  K <- nrow(u)
  D <- ncol(u)
  phi <- array(0, c(K, K, D + 1))
  for (d in seq_len(D)) {
    phi[, , d] <- rep(u[, d], each = K)  # phi[i, j, d] = u[j, d] (pure column)
  }
  phi[, , D + 1] <- tcrossprod(u)        # state similarity interaction
  if (!is.null(colnames(u))) {
    dimnames(phi) <- list(NULL, NULL, c(colnames(u), "sim"))
  }
  phi
}
