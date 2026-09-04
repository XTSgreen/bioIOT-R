# bioIOT single-cell interfaces (matrix / SingleCellExperiment / Seurat)

#' Run IOT on a single-cell object
#'
#' Aggregates cells into states (clusters), builds state-transition
#' features from the cell embedding, fits IOT weights, and returns the
#' state transition matrix with random-walk pseudotime. Generic with
#' methods for base matrices, \code{SingleCellExperiment} and
#' \code{Seurat} objects.
#'
#' @param object Cell embedding matrix (cells x dims) — for the matrix
#'   method — or a \code{SingleCellExperiment} / \code{Seurat} object.
#' @param state Factor/character of state (cluster) labels, one per cell
#'   (matrix method).
#' @param from,to Logical vectors or indices marking source / target
#'   timepoint cells (matrix method).
#' @param root Optional root state (name or index) for pseudotime.
#' @param ... Method-specific arguments.
#'
#' @return List with \code{fit} (\code{bioIOT_fit}), \code{Q} (transition
#'   matrix), \code{u} (state centroids), \code{a}, \code{b},
#'   \code{phi}, and \code{pseudotime} (if \code{root} given).
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_iot_states(K = 5, n_cells = 40, seed = 1)
#' res <- runIOT(sim$embedding, sim$cell_state,
#'               from = sim$cell_time == "t0", to = sim$cell_time == "t1",
#'               root = "S1", n_restart = 1, epochs = 100)
#' res$Q[1:3, 1:3]
#' @export
runIOT <- function(object, ...) {
  UseMethod("runIOT")
}

#' @rdname runIOT
#' @export
runIOT.default <- function(object, ...) {
  stop("runIOT() needs a cell embedding matrix, SingleCellExperiment ",
       "or Seurat object", call. = FALSE)
}

#' @rdname runIOT
#' @param n_dim Number of embedding dims used as state features.
#' @param T_obs Optional (K, K) observed row-conditional transition matrix
#'   (e.g. from lineage/clone data). When given, feature weights are fitted;
#'   otherwise the plan is solved with uniform feature weights.
#' @param mu,eps,lam,epochs,lr,n_restart,seed,two_stage Passed to
#'   \code{\link{fit_iot}} (only used when \code{T_obs} is given).
#' @importFrom stats runif
#' @export
runIOT.matrix <- function(object, state, from, to,
                          root = NULL, n_dim = NULL, T_obs = NULL,
                          mu = 0.5, eps = 1, lam = 0.05,
                          epochs = 300, lr = 0.01, n_restart = 4,
                          seed = 1, two_stage = TRUE, ...) {
  object <- as.matrix(object)
  state <- as.character(state)
  if (nrow(object) != length(state)) {
    stop("`state` must have one label per cell (row)", call. = FALSE)
  }
  if (is.logical(from)) from <- which(from)
  if (is.logical(to)) to <- which(to)
  states <- sort(unique(state[c(from, to)]))
  K <- length(states)
  D <- if (is.null(n_dim)) ncol(object) else min(n_dim, ncol(object))
  emb <- object[, seq_len(D), drop = FALSE]

  # state centroids in embedding space (all cells of the state)
  u <- t(vapply(states, function(s) colMeans(emb[state == s, , drop = FALSE]),
                numeric(D)))
  a_raw <- as.numeric(table(factor(state[from], levels = states)))
  b_raw <- as.numeric(table(factor(state[to], levels = states)))
  phi <- build_state_features(u)
  phi_z <- zscore_phi(phi)
  phi_use <- phi_z$phi_z

  if (!is.null(T_obs)) {
    T_obs <- as.matrix(T_obs)
    if (!identical(dim(T_obs), c(K, K))) {
      stop(sprintf("`T_obs` must be (%d, %d)", K, K), call. = FALSE)
    }
    fit <- fit_iot(list(phi_use), list(a_raw), list(b_raw), list(T_obs),
                   mu = mu, eps = eps, lam = lam, epochs = epochs, lr = lr,
                   n_restart = n_restart, seed = seed,
                   two_stage = two_stage, standardize = FALSE, ...)
    theta <- fit$theta
  } else {
    fit <- NULL
    theta <- rep(1, dim(phi_use)[3])
  }

  C <- make_cost(phi_use, theta)
  P <- soft_sinkhorn(C, a_raw, b_raw, mu = mu, eps = eps)
  a_norm <- a_raw / sum(a_raw)
  Q <- row_conditional(P, a_norm)
  rownames(Q) <- colnames(Q) <- states
  colnames(u) <- colnames(object)[seq_len(D)]
  rownames(u) <- states
  out <- list(fit = fit, Q = Q, u = u, a = a_norm,
              b = b_raw / sum(b_raw), phi = phi, phi_meta = phi_z$meta,
              states = states, theta = theta)
  if (!is.null(root)) {
    out$pseudotime <- pseudotime_from_transition(Q, root)
  }
  out
}

#' @rdname runIOT
#' @param state_col,colData_col State label column in colData (or a vector).
#' @param time_col,colData_time Timepoint column in colData (or a vector).
#' @param dimred Reduced dimension name (default "PCA").
#' @export
runIOT.SingleCellExperiment <- function(object, state_col, time_col,
                                        from, to, dimred = "PCA",
                                        root = NULL, ...) {
  stopifnot(requireNamespace("SingleCellExperiment", quietly = TRUE))
  emb <- as.matrix(SingleCellExperiment::reducedDim(object, dimred))
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  state <- if (state_col %in% names(cd)) cd[[state_col]] else state_col
  timev <- if (time_col %in% names(cd)) cd[[time_col]] else time_col
  res <- runIOT.matrix(emb, state = state,
                       from = timev == from, to = timev == to,
                       root = root, ...)
  res
}

#' @rdname runIOT
#' @param group.by Seurat ident / metadata column for states.
#' @param split.by Metadata column holding the timepoint values
#'   \code{from}/\code{to}.
#' @param reduction Reduction name (default "pca").
#' @export
runIOT.Seurat <- function(object, group.by = NULL, split.by = NULL,
                          from = NULL, to = NULL, reduction = "pca",
                          root = NULL, ...) {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  emb <- as.matrix(Seurat::Embeddings(object, reduction = reduction))
  meta <- object[[]]
  state <- if (!is.null(group.by)) {
    if (group.by == "ident") as.character(Seurat::Idents(object)) else meta[[group.by]]
  } else {
    as.character(Seurat::Idents(object))
  }
  timev <- if (!is.null(split.by)) meta[[split.by]] else NULL
  if (is.null(timev)) {
    stop("`split.by` with `from`/`to` values is required", call. = FALSE)
  }
  res <- runIOT.matrix(emb, state = state,
                       from = timev == from, to = timev == to,
                       root = root, ...)
  res
}
