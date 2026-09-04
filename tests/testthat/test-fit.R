test_that("fit_iot recovers a single pure-column feature", {
  K <- 6
  set.seed(6)
  a <- rep(1, K) / K
  b <- rep(1, K) / K
  u <- runif(K, 0.3, 3)
  phi <- array(rep(u, each = K), c(K, K, 1))
  phi_z <- zscore_phi(phi)$phi_z
  T_obs <- row_conditional(
    soft_sinkhorn(make_cost(phi_z, 2.0), a, b, mu = 0.5), a)
  fit <- fit_iot(phi_z, a, b, T_obs, mu = 0.5, n_restart = 1,
                 seed = 1, epochs = 200)
  expect_true(is.finite(fit$theta))
  expect_gt(fit$theta, 0.5)
  expect_lt(fit$theta, 4)
  expect_true(is.numeric(fit$loss))
})

test_that("fit_iot: class, multi-scenario list input, print/summary", {
  K <- 4
  F <- 2
  set.seed(7)
  a <- rep(1, K) / K
  b <- rep(1, K) / K
  phi1 <- array(rnorm(K * K * F), c(K, K, F))
  phi2 <- array(rnorm(K * K * F), c(K, K, F))
  T1 <- matrix(runif(K * K), K); T1 <- sweep(T1, 1, rowSums(T1), "/")
  T2 <- matrix(runif(K * K), K); T2 <- sweep(T2, 1, rowSums(T2), "/")
  fit <- fit_iot(list(phi1, phi2), list(a, a), list(b, b), list(T1, T2),
                 n_restart = 1, epochs = 5, seed = 0)
  expect_s3_class(fit, "bioIOT_fit")
  expect_equal(length(fit$scenarios), 2)
  expect_equal(length(fit$theta), F)
  expect_true(is.logical(fit$support))
  expect_equal(length(fit$restart_losses), 1)
  expect_output(print(fit), "bioIOT fit")
  expect_output(print(summary(fit)), "bioIOT fit")
})

test_that("fit_iot: input validation", {
  K <- 3
  expect_error(fit_iot(array(1, c(K, K, 2)), rep(1, K), rep(1, K),
                       matrix(1, K, K)), NA)
  expect_error(fit_iot(array(1, c(K, K, 2)), rep(1, K + 1), rep(1, K),
                       matrix(1, K, K)), "length")
  expect_error(fit_iot(list(array(1, c(K, K, 2)), array(1, c(4, 4, 2))),
                       rep(1, K), rep(1, K), matrix(1, K, K)), "scenarios")
})

test_that("fit_iot: scenarios may differ in K (empty states at some timepoints)", {
  set.seed(11)
  K1 <- 5
  a1 <- rep(1, K1) / K1
  phi1 <- array(rnorm(K1 * K1 * 2), c(K1, K1, 2))
  T1 <- matrix(runif(K1 * K1), K1); T1 <- sweep(T1, 1, rowSums(T1), "/")
  K2 <- 6
  a2 <- rep(1, K2) / K2
  phi2 <- array(rnorm(K2 * K2 * 2), c(K2, K2, 2))
  T2 <- matrix(runif(K2 * K2), K2); T2 <- sweep(T2, 1, rowSums(T2), "/")
  fit <- fit_iot(list(phi1, phi2), list(a1, a2), list(a1, a2), list(T1, T2),
                 n_restart = 1, epochs = 10, seed = 0)
  expect_equal(fit$K, c(5, 6))
  expect_equal(length(fit$theta), 2)
  Q2 <- transition_matrix(fit, which = 2)
  expect_equal(dim(Q2), c(6, 6))
  expect_error(fit_iot(list(phi1, array(1, c(4, 4, 3))), list(a1, a2),
                       list(a1, a2), list(T1, T2)), "K, K, 2")
})
