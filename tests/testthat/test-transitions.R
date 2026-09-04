test_that("transition_matrix: row-stochastic, from stored scenario or override", {
  set.seed(1)
  sim <- simulate_iot_states(K = 5, seed = 1)
  fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T_true,
                 n_restart = 1, epochs = 60, seed = 1)
  Q <- transition_matrix(fit)
  expect_equal(dim(Q), c(5, 5))
  expect_equal(unname(rowSums(Q)), rep(1, 5), tolerance = 1e-8)
  # override path
  Q2 <- transition_matrix(fit, phi = sim$phi, a = sim$a, b = sim$b)
  expect_equal(dim(Q2), c(5, 5))
  expect_error(transition_matrix(list()), "bioIOT_fit")
})

test_that("pseudotime_from_transition: root zero, others positive", {
  Q <- matrix(c(0.5, 0.5, 0, 0,
                0.2, 0.8, 0, 0,
                0.1, 0.1, 0.4, 0.4,
                0, 0.1, 0.1, 0.8), 4, 4, byrow = TRUE)
  pt <- pseudotime_from_transition(Q, root = 1)
  expect_equal(unname(pt[1]), 0)
  expect_true(all(pt[-1] > 0))
  expect_gt(pt[2], pt[1])       # S2 directly feeds S1
  expect_error(pseudotime_from_transition(Q, root = "nope"), "root")
  pt2 <- pseudotime_from_transition(Q, root = "S1")  # unnamed -> S-labels
  expect_equal(unname(pt2[1]), 0)
})

test_that("build_state_features: pure-column + similarity blocks", {
  u <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2, dimnames = list(NULL, c("A", "B")))
  phi <- build_state_features(u)
  expect_equal(dim(phi), c(3, 3, 3))
  expect_equal(dimnames(phi)[[3]], c("A", "B", "sim"))
  expect_equal(as.vector(phi[, , 1]), rep(u[, 1], each = 3))  # phi[i,j] = u[j]
  expect_equal(as.vector(phi[, , 3]), as.vector(tcrossprod(u)))
  expect_error(build_state_features(matrix(c(NA, 1), 2, 1)), "non-finite")
})

test_that("pseudotime_from_transition: zero-outflow terminal states get late times", {
  Q <- matrix(c(0.5, 0.5, 0, 0,
                0.2, 0.8, 0, 0,
                0.1, 0.1, 0.4, 0.4,
                0, 0, 0, 0), 4, 4, byrow = TRUE)  # S4 terminal (zero row)
  pt <- pseudotime_from_transition(Q, root = 1)
  expect_equal(unname(pt[1]), 0)
  expect_true(all(is.finite(pt)))
  expect_gt(pt[4], min(pt[-1]))          # terminal state is late, not transient
})
