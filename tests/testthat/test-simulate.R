test_that("simulate_iot_states: reproducible, sane shapes", {
  sim1 <- simulate_iot_states(K = 5, n_cells = 20, seed = 1)
  sim2 <- simulate_iot_states(K = 5, n_cells = 20, seed = 1)
  expect_equal(sim1$P_true, sim2$P_true, tolerance = 1e-12)
  expect_equal(dim(sim1$phi), c(5, 5, 3))
  expect_equal(dim(sim1$cell_embedding), c(200, 2))
  expect_equal(length(sim1$cell_state), 200)
  expect_equal(unname(rowSums(sim1$T_true)), rep(1, 5), tolerance = 1e-8)
})

test_that("demo data ships and loads", {
  skip_on_cran()
  data(demo_iot_states, package = "bioIOT", envir = environment())
  expect_equal(demo_iot_states$K, 6)
  expect_true(is.matrix(demo_iot_states$P_true))
})
