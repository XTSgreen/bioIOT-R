# bioIOT simulated demo data generator

#' Simulate single-cell state-transition data for bioIOT
#'
#' Generates a reproducible synthetic single-cell-like dataset: states with
#' feature centroids, source/target timepoint masses, a true IOT plan and
#' transition matrix, plus cell-level metadata for the \code{\link{runIOT}}
#' interface. The bundled \code{\link{demo_iot_states}} dataset is produced
#' by this function with its default seed.
#'
#' @param K Number of states.
#' @param F_emb Number of state feature dimensions (>= 2 for plotting).
#' @param n_cells Cells per state per timepoint.
#' @param seed Random seed.
#' @param theta_true True feature weights (length F_emb + 1).
#' @param mu,eps Solver hyperparameters used to build the true plan.
#'
#' @return List with \code{u} (K x F_emb state centroids), \code{phi}
#'   ((K, K, F_emb + 1) features), \code{a}, \code{b} (state masses),
#'   \code{P_true}, \code{T_true} (row-conditional), \code{theta_true},
#'   \code{embedding} (K x 2 state coordinates for plotting),
#'   \code{cell_state} / \code{cell_time} / \code{cell_embedding}
#'   (cell-level metadata for \code{runIOT}).
#'
#' @examples
#' sim <- simulate_iot_states(K = 5, seed = 1)
#' names(sim)
#' @export
simulate_iot_states <- function(K = 6, F_emb = 2, n_cells = 50, seed = 1,
                                theta_true = c(0.9, -0.7, 1.1),
                                mu = 0.5, eps = 1) {
  set.seed(seed)
  u <- matrix(runif(K * F_emb, 0.3, 3), K, F_emb)
  phi <- build_state_features(u)
  phi_z <- zscore_phi(phi)
  if (length(theta_true) != dim(phi)[3]) {
    stop(sprintf("`theta_true` must have length %d (F_emb + 1)", dim(phi)[3]),
         call. = FALSE)
  }
  a_raw <- rep(1, K)
  b_raw <- 1 + 0.5 * (seq_len(K) %% 2)  # non-uniform target composition
  a <- a_raw / sum(a_raw)
  b <- b_raw / sum(b_raw)
  C <- make_cost(phi_z$phi_z, theta_true)
  P <- soft_sinkhorn(C, a, b, mu = mu, eps = eps)
  T_true <- row_conditional(P, a)

  # 2-D embedding for plots: first two feature dims with slight jitter
  emb <- u[, seq_len(min(2, F_emb)), drop = FALSE] +
    matrix(runif(K * min(2, F_emb), -0.15, 0.15), K)

  # cell-level metadata: two timepoints, cells jittered around centroids
  D <- F_emb
  cell_emb <- do.call(rbind, lapply(seq_len(K), function(k) {
    m1 <- c(u[k, ], rep(0, max(0, D - F_emb)))[seq_len(D)]
    sapply(seq_len(D), function(d) rnorm(n_cells * 2, m1[d], 0.25))
  }))
  cell_state <- factor(rep(sprintf("S%d", seq_len(K)), each = n_cells * 2))
  cell_time <- factor(rep(rep(c("t0", "t1"), each = n_cells), K))

  list(u = u, phi = phi, a = a / sum(a), b = b / sum(b),
       P_true = P, T_true = T_true, theta_true = theta_true,
       embedding = emb, K = K, F_emb = F_emb,
       cell_embedding = cell_emb, cell_state = cell_state,
       cell_time = cell_time)
}
