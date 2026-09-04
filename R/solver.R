# bioIOT core solver
#
# 半松弛 OT 求解器的 R 实现（与 Python 包/论文核心数值同源）：
#   min_P <C,P> - eps*H(P) + mu*KL(col(P)||b),  s.t. P1 = a
# 前向：Anderson 加速阻尼不动点迭代；反向：不动点隐式微分。

#' Coerce to a normalized probability vector
#' @noRd
.prob_vector <- function(x, size, name) {
  x <- as.numeric(x)
  if (length(x) != size) {
    stop(sprintf("`%s` must have length %d, got %d", name, size, length(x)),
         call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("`%s` must be finite and non-negative", name), call. = FALSE)
  }
  s <- sum(x)
  if (s <= 0) {
    stop(sprintf("`%s` must have positive total mass", name), call. = FALSE)
  }
  x / s
}

#' Rebuild the plan P from a converged column marginal (col = P^T 1)
#' @noRd
.plan_from_col <- function(col, C, a, b, mu, eps) {
  K <- nrow(C)
  logits <- -C / eps +
    (mu / eps) * matrix(log(b / (col + 1e-300)), K, K, byrow = TRUE)
  logits <- sweep(logits, 1, apply(logits, 1, max), "-")
  e <- exp(logits)
  s <- e / rowSums(e)
  a * s  # column-major recycling: P_ij = a_i * s_ij
}

#' One step of the column-marginal fixed-point map
#' @noRd
.fcol <- function(col, C, a, b, mu, eps) {
  K <- nrow(C)
  logits <- -C / eps +
    (mu / eps) * matrix(log(b / (col + 1e-300)), K, K, byrow = TRUE)
  logits <- sweep(logits, 1, apply(logits, 1, max), "-")
  e <- exp(logits)
  s <- e / rowSums(e)
  as.vector(crossprod(s, a))
}

#' Anderson-accelerated damped fixed-point solve for the column marginal
#' @noRd
.anderson_col <- function(col0, C, a, b, mu, eps, iters, damp,
                          m = 6L, tol = 1e-12) {
  col <- col0
  xs <- list()
  gs <- list()
  best_col <- col
  best_r <- Inf
  for (it in seq_len(iters)) {
    fx <- .fcol(col, C, a, b, mu, eps)
    g <- fx - col
    r <- max(abs(g))
    if (r < best_r) {
      best_r <- r
      best_col <- fx
    }
    if (r < tol) {
      return(list(col = fx, P = .plan_from_col(fx, C, a, b, mu, eps)))
    }
    xs[[length(xs) + 1L]] <- col
    gs[[length(gs) + 1L]] <- g
    if (length(xs) > m) {
      xs[[1L]] <- NULL
      gs[[1L]] <- NULL
    }
    col_dmp <- damp * col + (1 - damp) * fx
    if (length(xs) >= 2L) {
      col_and <- tryCatch(
        {
          gm <- gs[[length(gs)]]
          xm <- xs[[length(xs)]]
          G <- vapply(gs[-length(gs)], function(gg) gg - gm, numeric(length(gm)))
          if (is.null(dim(G))) G <- matrix(G, nrow = length(gm))
          coef <- -qr.solve(G, gm)
          X <- vapply(seq_along(xs[-length(xs)]),
                      function(i) (xs[[i]] - xm) + (gs[[i]] - gm),
                      numeric(length(xm)))
          if (is.null(dim(X))) X <- matrix(X, nrow = length(xm))
          cand <- xm + gm + as.vector(X %*% coef)
          if (all(is.finite(cand))) {
            cand <- pmax(cand, 1e-12)
            cand <- cand / sum(cand)
            fx_and <- .fcol(cand, C, a, b, mu, eps)
            r_and <- max(abs(fx_and - cand))
            r_dmp <- max(abs(.fcol(col_dmp, C, a, b, mu, eps) - col_dmp))
            if (r_and < r_dmp) cand else col_dmp
          } else {
            col_dmp
          }
        },
        error = function(e) col_dmp,
        warning = function(w) col_dmp
      )
      col <- col_and
    } else {
      col <- col_dmp
    }
  }
  if (best_r > 1e-4) {
    warning(sprintf("[bioIOT] fixed-point residual %.2e after %d iterations",
                    best_r, iters), call. = FALSE)
  }
  list(col = best_col, P = .plan_from_col(best_col, C, a, b, mu, eps))
}

#' Solve the semi-relaxed OT problem
#'
#' Solves
#' \deqn{\min_P <C,P> - \eps H(P) + \mu KL(col(P) || b)\ \ \mathrm{s.t.}\ P\,\mathbf{1}=a}{}
#' where the source-side row marginal is hard and the column marginal is
#' KL-anchored to \code{b} with strength \code{mu}.
#'
#' @param C (K, K) cost matrix.
#' @param a,b (K,) source / target masses; any positive scale is accepted
#'   (auto-normalized). Zero-mass states are allowed in \code{a}.
#' @param mu Soft-marginal strength; large \code{mu} approaches hard-marginal
#'   OT, \code{mu} = 0 recovers a plain row-softmax.
#' @param eps Entropic regularization (> 0).
#' @param iters Fixed-point iteration cap.
#' @param damp Damping factor of the base iteration.
#' @param tol Fixed-point residual tolerance.
#'
#' @return (K, K) transport plan with row sums equal to the normalized
#'   \code{a}.
#' @examples
#' K <- 5
#' set.seed(1)
#' P <- soft_sinkhorn(matrix(rnorm(25), 5), rep(1, 5), rep(1, 5))
#' rowSums(P)
#' @export
soft_sinkhorn <- function(C, a, b, mu = 0.5, eps = 1,
                          iters = 1000L, damp = 0.5, tol = 1e-12) {
  C <- as.matrix(C)
  if (!is.numeric(C)) {
    C <- matrix(as.numeric(C), nrow(C), ncol(C))
  }
  if (any(!is.finite(C))) stop("`C` contains non-finite values", call. = FALSE)
  K <- nrow(C)
  if (ncol(C) != K) {
    stop(sprintf("`C` must be square, got %d x %d", K, ncol(C)), call. = FALSE)
  }
  a <- .prob_vector(a, K, "a")
  b <- .prob_vector(b, K, "b")
  if (eps <= 0) stop("`eps` must be > 0", call. = FALSE)
  if (mu < 0) stop("`mu` must be >= 0", call. = FALSE)
  .anderson_col(b, C, a, b, mu, eps, iters, damp, tol = tol)$P
}

#' Row-conditional transition matrix
#'
#' \eqn{Q(i,\cdot) = P(i,\cdot) / a_i}: the state-to-state transition
#' probabilities conditional on starting in state \eqn{i}.
#'
#' @param P (K, K) transport plan.
#' @param a (K,) source masses (row sums of \code{P}; zeros allowed).
#' @param eps_t Numerical floor.
#'
#' @return (K, K) matrix; rows with \eqn{a_i = 0} are all-zero.
#' @examples
#' P <- matrix(1, 3, 3) / 9
#' row_conditional(P, c(1, 0, 2) / 3)
#' @export
row_conditional <- function(P, a, eps_t = 1e-300) {
  P <- as.matrix(P)
  a <- as.numeric(a)
  if (length(a) != nrow(P)) {
    stop("`a` length must equal nrow(P)", call. = FALSE)
  }
  if (any(!is.finite(a)) || any(a < 0)) {
    stop("`a` must be finite and non-negative", call. = FALSE)
  }
  Q <- P / (a + eps_t)
  zero <- a <= eps_t
  if (any(zero)) Q[zero, ] <- 0  # zero-mass source rows have no transition
  Q
}
