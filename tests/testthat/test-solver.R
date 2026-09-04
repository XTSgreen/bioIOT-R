test_that("soft_sinkhorn: row marginal exact, plan finite non-negative", {
  K <- 8
  set.seed(0)
  C <- matrix(rnorm(K * K), K)
  a <- rep(1, K) / K
  b <- rep(1, K) / K
  P <- soft_sinkhorn(C, a, b, mu = 0.5)
  expect_true(all(is.finite(P)))
  expect_gte(min(P), 0)
  expect_equal(rowSums(P), a, tolerance = 1e-8)
})

test_that("soft_sinkhorn: unnormalized masses accepted (scale invariant)", {
  K <- 5
  set.seed(1)
  C <- matrix(rnorm(K * K), K)
  a_raw <- rep(37, K)
  b_raw <- (1:K) * 11
  a <- a_raw / sum(a_raw)
  b <- b_raw / sum(b_raw)
  P1 <- soft_sinkhorn(C, a, b, mu = 0.5)
  P2 <- soft_sinkhorn(C, a_raw, b_raw, mu = 0.5)
  expect_equal(P1, P2, tolerance = 1e-10)
  expect_equal(rowSums(P2), a, tolerance = 1e-8)
})

test_that("soft_sinkhorn: larger mu anchors col(P) closer to b", {
  K <- 6
  set.seed(2)
  C <- matrix(rnorm(K * K), K)
  a <- rep(1, K) / K
  b <- exp(seq(0, 0.7, length.out = K)); b <- b / sum(b)
  dev <- vapply(c(0.05, 5), function(mu) {
    P <- soft_sinkhorn(C, a, b, mu = mu, iters = 5000)
    max(abs(colSums(P) - b))
  }, numeric(1))
  expect_lt(dev[2], dev[0 + 1])
})

test_that("soft_sinkhorn: input validation", {
  K <- 4
  expect_error(soft_sinkhorn(matrix(0, K, K + 1), rep(1, K), rep(1, K)), "square")
  expect_error(soft_sinkhorn(matrix(1, K, K), rep(-1, K), rep(1, K)), "non-negative")
  expect_error(soft_sinkhorn(matrix(1, K, K), rep(0, K), rep(1, K)), "positive total mass")
  expect_error(soft_sinkhorn(matrix(1, K, K), rep(1, K), rep(1, K), eps = 0), "eps")
  expect_error(soft_sinkhorn(matrix(NA, K, K), rep(1, K), rep(1, K)), "non-finite")
})

test_that("row_conditional: zero-mass rows are zero", {
  P <- matrix(1, 3, 3)
  a <- c(1, 0, 2) / 3
  Q <- row_conditional(P, a)
  expect_equal(sum(Q[2, ]), 0)
  expect_equal(Q[3, ], rep(1.5, 3), tolerance = 1e-9)
})
