# bioIOT fitting: implicit-gradient Adam + two-stage debias + multi-restart
#
# 与论文/Python 包 fit_uot 同构：
#   stage 1: l1 正则 Adam（隐式微分梯度）
#   stage 2: 支撑集上去偏重拟合（Bühlmann 式）
#   多重启取最终无正则损失最低者。

#' Gradient of the plan w.r.t. the cost via the implicit function theorem
#' @noRd
.dLdC <- function(P, col, a, b, mu, eps, damp, gP) {
  K <- nrow(P)
  s <- P / (a + 1e-300)                       # row softmax values s_ij
  Tmat <- crossprod(P, s)                     # T[j,k] = sum_i P_ij s_ik
  invcol <- 1 / (col + 1e-300)
  A <- diag(damp, K) -
    (1 - damp) * (mu / eps) * sweep(diag(col) - Tmat, 2, invcol, `*`)
  rowsum_g <- rowSums(gP * s)
  direct <- -(a / eps) * s * (gP - rowsum_g)
  term1 <- colSums(gP * P)
  rowsum_gP <- rowSums(gP * P)
  term2 <- as.vector(crossprod(s, rowsum_gP))
  w <- -(mu / eps) * invcol * (term1 - term2)
  lhs <- diag(K) - t(A) + diag(1e-12, K)
  adj <- as.vector(solve(lhs, w))
  adj_rowsum <- as.vector(s %*% adj)
  dmat <- matrix(adj, K, K, byrow = TRUE) - adj_rowsum
  corr <- -(1 - damp) * (a / eps) * s * dmat
  direct + corr
}

#' Gradient of the row-CE loss w.r.t. the plan
#' @noRd
.grad_P_ce <- function(T, P, a, eps_t = 1e-12) {
  Q <- P / (a + 1e-300)
  mask <- which(a > eps_t)
  G <- matrix(0, nrow(P), ncol(P))
  G[mask, ] <- -T[mask, , drop = FALSE] /
    ((Q[mask, , drop = FALSE] + eps_t) * (a[mask] + 1e-300))
  G
}

#' Gradient of the row-CE loss w.r.t. theta (all scenarios summed)
#' @noRd
.grad_theta <- function(scen, theta, mu, eps, iters, damp) {
  loss <- 0
  grad <- rep(0, length(theta))
  for (sc in scen) {
    C <- make_cost(sc$phi, theta)
    res <- .anderson_col(sc$b, C, sc$a, sc$b, mu, eps, iters, damp)
    gP <- .grad_P_ce(sc$T, res$P, sc$a)
    dLdC <- .dLdC(res$P, res$col, sc$a, sc$b, mu, eps, damp, gP)
    loss <- loss + row_ce_loss(sc$T, res$P, sc$a)
    # C = -einsum(phi, theta)  =>  dL/dtheta_k = -sum_ij dLdC_ij * phi_ijk
    grad <- grad - vapply(seq_len(dim(sc$phi)[3]),
                          function(k) sum(dLdC * sc$phi[, , k]), numeric(1))
  }
  list(loss = loss, grad = grad)
}

#' Adam loop (paper parity: warm-up then l1, grad clip, finite restart)
#' @noRd
.adam_loop <- function(theta0, loss_grad_fn, epochs, lr, lam, warm_frac = 0.4) {
  theta <- theta0
  m <- rep(0, length(theta))
  v <- rep(0, length(theta))
  b1 <- 0.9
  b2 <- 0.999
  e <- 1e-8
  warm <- floor(epochs * warm_frac)
  for (it in seq_len(epochs)) {
    lg <- loss_grad_fn(theta)
    g <- lg$grad
    reg <- if (it > warm) lam * sign(theta) else rep(0, length(theta))
    g <- g + reg
    gn <- sqrt(sum(g^2))  # global L2 clip to norm 1 (paper parity)
    if (gn > 1) g <- g / gn
    m <- b1 * m + (1 - b1) * g
    v <- b2 * v + (1 - b2) * g^2
    mhat <- m / (1 - b1^it)
    vhat <- v / (1 - b2^it)
    theta <- pmin(pmax(theta - lr * mhat / (sqrt(vhat) + e), -10), 10)
    if (!all(is.finite(theta))) theta <- rnorm(length(theta), 0, 0.05)
  }
  theta
}

#' Validate and normalize a set of scenarios
#'
#' Scenarios may differ in the number of states K (e.g. empty clusters at
#' some timepoints); the number of features F must be shared.
#' @noRd
.check_scenarios <- function(phil, a_l, b_l, T_l) {
  if (!(length(phil) == length(a_l) && length(a_l) == length(b_l) &&
        length(b_l) == length(T_l))) {
    stop("phil, a_l, b_l, T_l must have the same number of scenarios",
         call. = FALSE)
  }
  dims0 <- dim(phil[[1]])
  if (length(dims0) != 3L) stop("each `phi` must be (K, K, F)", call. = FALSE)
  F0 <- dims0[3]
  scen <- vector("list", length(phil))
  for (i in seq_along(phil)) {
    d_i <- dim(phil[[i]])
    if (is.null(d_i) || length(d_i) != 3L || d_i[3] != F0) {
      stop(sprintf("scenario %d: `phi` must be (K, K, %d)", i, F0),
           call. = FALSE)
    }
    if (d_i[1] != d_i[2]) {
      stop(sprintf("scenario %d: `phi` must be square (K, K, F)", i),
           call. = FALSE)
    }
    phi <- array(as.numeric(phil[[i]]), d_i)
    if (any(!is.finite(phi))) {
      stop(sprintf("scenario %d: `phi` contains non-finite values", i),
           call. = FALSE)
    }
    K <- d_i[1]
    a <- .prob_vector(a_l[[i]], K, "a")
    b <- .prob_vector(b_l[[i]], K, "b")
    T <- as.matrix(T_l[[i]])
    if (!identical(dim(T), c(K, K))) {
      stop(sprintf("scenario %d: `T` must be (%d, %d)", i, K, K), call. = FALSE)
    }
    if (any(!is.finite(T)) || any(T < 0)) {
      stop(sprintf("scenario %d: `T` must be finite and non-negative", i),
           call. = FALSE)
    }
    scen[[i]] <- list(phi = phi, a = a, b = b, T = T)
  }
  scen
}

#' Fit inverse optimal transport feature weights
#'
#' Fits feature weights \eqn{\theta} so that the soft-marginal OT plan
#' induced by \eqn{C = -\mathrm{einsum}(\phi, \theta)} reproduces the
#' observed row-conditional transitions \eqn{T}. Two-stage (l1 selection
#' then debias refit) with multi-restart; gradients are exact implicit
#' differentiations of the fixed point.
#'
#' @param phi (K, K, F) feature array, or a list of arrays (scenarios).
#' @param a,b (K,) masses, or lists (auto-normalized).
#' @param T (K, K) observed row-conditional transitions, or a list.
#' @param mu Soft-marginal strength (default 0.5, paper working point).
#' @param lam l1 strength in stage 1 (default 0.05).
#' @param eps Entropic regularization (default 1).
#' @param epochs Adam steps per stage (default 300).
#' @param lr Adam learning rate (default 0.01).
#' @param n_restart Random restarts (default 4).
#' @param seed Base seed for restart inits.
#' @param two_stage Debias refit on the selected support (default TRUE).
#' @param iters,damp Forward solver controls.
#' @param standardize Z-score each scenario's features before fitting
#'   (default TRUE, paper pipeline).
#' @param verbose Print progress.
#'
#' @return Object of class \code{bioIOT_fit} with components \code{theta}
#'   (debiased weights), \code{theta1} (stage-1 weights), \code{support}
#'   (logical mask), \code{loss} (final row-CE), \code{restart_losses},
#'   \code{scenarios} (stored standardized scenarios), \code{meta}
#'   (standardization transforms) and \code{hyper}.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_iot_states(K = 5, seed = 1)
#' fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
#' fit$theta
#' @export
fit_iot <- function(phi, a, b, T,
                    mu = 0.5, lam = 0.05, eps = 1,
                    epochs = 300, lr = 0.01,
                    n_restart = 4, seed = 1,
                    two_stage = TRUE, iters = 1000L, damp = 0.5,
                    standardize = TRUE, verbose = FALSE) {
  phil <- if (is.list(phi)) phi else list(phi)
  a_l <- if (is.list(a)) a else list(a)
  b_l <- if (is.list(b)) b else list(b)
  T_l <- if (is.list(T)) T else list(T)
  scen <- .check_scenarios(phil, a_l, b_l, T_l)
  F <- dim(scen[[1]]$phi)[3]

  meta <- vector("list", length(scen))
  if (standardize) {
    for (i in seq_along(scen)) {
      z <- zscore_phi(scen[[i]]$phi)
      scen[[i]]$phi <- z$phi_z
      meta[[i]] <- z$meta
    }
  }

  ce_grad <- function(theta) .grad_theta(scen, theta, mu, eps, iters, damp)

  restart_losses <- numeric(n_restart)
  best <- NULL
  for (r in seq_len(n_restart)) {
    set.seed(seed * 100 + r)
    theta1 <- .adam_loop(rnorm(F, 0, 0.05), ce_grad, epochs, lr, lam)
    support <- abs(theta1) > 1e-3
    theta_final <- theta1
    if (two_stage && any(support)) {
      theta2 <- theta1 * support
      theta_final <- .adam_loop(theta2, ce_grad,
                                max(1L, floor(epochs * 0.75)), lr * 2, 0)
    }
    loss_f <- ce_grad(theta_final)$loss
    restart_losses[r] <- loss_f
    if (is.null(best) || loss_f < best$loss) {
      best <- list(theta1 = theta1, support = support,
                   theta = theta_final, loss = loss_f)
    }
    if (verbose) {
      message(sprintf("[bioIOT] restart %d/%d loss=%.4f", r, n_restart, loss_f))
    }
  }

  out <- list(
    theta = best$theta,
    theta1 = best$theta1,
    support = best$support,
    loss = best$loss,
    restart_losses = restart_losses,
    scenarios = scen,
    meta = meta,
    hyper = list(mu = mu, lam = lam, eps = eps, epochs = epochs, lr = lr,
                 n_restart = n_restart, seed = seed, two_stage = two_stage,
                 iters = iters, damp = damp, standardize = standardize),
    K = vapply(scen, function(s) dim(s$phi)[1], integer(1)),
    F = F,
    call = match.call()
  )
  class(out) <- "bioIOT_fit"
  out
}

#' Print a bioIOT fit
#' @export
#' @param x A \code{bioIOT_fit} object.
#' @param ... Unused.
print.bioIOT_fit <- function(x, ...) {
  ks <- unique(x$K)
  ktxt <- if (length(ks) == 1L) sprintf("K=%d", ks) else
    sprintf("K=%s", paste(range(ks), collapse = "-"))
  cat(sprintf("bioIOT fit: %s states, F=%d features, %d scenario(s)\n",
              ktxt, x$F, length(x$scenarios)))
  cat(sprintf("  final row-CE loss: %.4f (%d restart(s))\n",
              x$loss, length(x$restart_losses)))
  cat(sprintf("  support: %d / %d features selected\n", sum(x$support), x$F))
  invisible(x)
}

#' Summarize a bioIOT fit
#' @export
#' @param object A \code{bioIOT_fit} object.
summary.bioIOT_fit <- function(object, ...) {
  out <- data.frame(
    feature = seq_len(object$F),
    theta = round(object$theta, 4),
    theta1 = round(object$theta1, 4),
    support = object$support,
    row.names = NULL
  )
  cat(sprintf("bioIOT fit (loss %.4f)\n", object$loss))
  print(out)
  invisible(out)
}
