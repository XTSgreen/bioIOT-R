test_that("implicit gradient matches central finite differences", {
  K <- 5
  F <- 2
  set.seed(3)
  a <- rep(1, K) / K
  b <- rep(1, K) / K
  phi <- array(rnorm(K * K * F), c(K, K, F))
  theta_star <- c(0.9, -1.3)
  P0 <- soft_sinkhorn(make_cost(phi, theta_star), a, b, mu = 0.5)
  T_obs <- row_conditional(P0, a)
  th0 <- c(0.4, -0.5)

  L <- function(th) row_ce_loss(T_obs, soft_sinkhorn(make_cost(phi, th), a, b, mu = 0.5), a)

  g_analytic <- bioIOT:::.grad_theta(
    list(list(phi = phi, a = a, b = b, T = T_obs)),
    th0, mu = 0.5, eps = 1, iters = 1000L, damp = 0.5
  )$grad

  h <- 1e-4
  g_fd <- vapply(seq_len(F), function(k) {
    e <- rep(0, F); e[k] <- h
    (L(th0 + e) - L(th0 - e)) / (2 * h)
  }, numeric(1))

  expect_equal(g_analytic, g_fd, tolerance = 1e-3)
})
