test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  Q <- matrix(c(0.7, 0.3, 0, 0.2, 0.6, 0.2, 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
  emb <- matrix(c(0, 0, 1, 1, 2, 0), 3, 2, byrow = TRUE)
  expect_true(ggplot2::is_ggplot(plot_transition_heatmap(Q)))
  expect_true(ggplot2::is_ggplot(plot_transition_flow(Q, emb)))
  expect_true(ggplot2::is_ggplot(plot_transition_flow(Q, emb, threshold = 0.9)))

  set.seed(1)
  sim <- simulate_iot_states(K = 5, seed = 1)
  fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T_true,
                 n_restart = 1, epochs = 60, seed = 1)
  expect_true(ggplot2::is_ggplot(plot_theta(fit)))
  expect_true(ggplot2::is_ggplot(plot_theta(sim$theta_true)))

  expr <- matrix(rnorm(200 * 12), nrow = 200)
  rownames(expr) <- paste0("G", 1:200)
  colnames(expr) <- paste0("S", 1:12)
  rownames(expr)[1] <- "VIM"; rownames(expr)[2] <- "CDH2"
  pw <- score_pathways(expr)
  expect_true(ggplot2::is_ggplot(plot_pathway_trend(pw, time = sort(runif(12)))))

  expect_error(plot_transition_flow(Q, emb[, 1, drop = FALSE]), "2 columns")
})
