test_that("runIOT.matrix: unsupervised plan + pseudotime", {
  set.seed(1)
  sim <- simulate_iot_states(K = 5, n_cells = 40, seed = 1)
  res <- runIOT(sim$cell_embedding, sim$cell_state,
                from = sim$cell_time == "t0", to = sim$cell_time == "t1",
                root = "S1")
  expect_null(res$fit)
  expect_equal(dim(res$Q), c(5, 5))
  expect_equal(unname(rowSums(res$Q)), rep(1, 5), tolerance = 1e-8)
  expect_equal(unname(res$pseudotime["S1"]), 0)
  expect_true(all(res$pseudotime[-1] > 0))
})

test_that("runIOT.matrix: supervised mode with T_obs fits theta", {
  set.seed(2)
  sim <- simulate_iot_states(K = 5, n_cells = 40, seed = 2)
  res <- runIOT(sim$cell_embedding, sim$cell_state,
                from = sim$cell_time == "t0", to = sim$cell_time == "t1",
                T_obs = sim$T_true, n_restart = 1, epochs = 80, seed = 1)
  expect_s3_class(res$fit, "bioIOT_fit")
  expect_true(all(is.finite(res$theta)))
})

test_that("runIOT.SingleCellExperiment works on a tiny object", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  set.seed(3)
  sim <- simulate_iot_states(K = 4, n_cells = 15, seed = 3)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(rnorm(20 * 120), 20, 120)),
    colData = S4Vectors::DataFrame(state = sim$cell_state,
                                   time = sim$cell_time)
  )
  SingleCellExperiment::reducedDim(sce, "PCA") <- as.matrix(sim$cell_embedding)
  res <- runIOT(sce, state_col = "state", time_col = "time",
                from = "t0", to = "t1", dimred = "PCA", root = "S1",
                n_dim = 2)
  expect_equal(dim(res$Q), c(4, 4))
  expect_equal(unname(rowSums(res$Q)), rep(1, 4), tolerance = 1e-8)
})

test_that("runIOT.Seurat works on a tiny object", {
  skip_if_not_installed("Seurat")
  set.seed(4)
  sim <- simulate_iot_states(K = 4, n_cells = 15, seed = 4)
  counts <- matrix(rpois(20 * 120, 5), 20, 120,
                   dimnames = list(paste0("g", 1:20), paste0("c", 1:120)))
  obj <- Seurat::CreateSeuratObject(counts = counts,
                                    meta.data = data.frame(
                                      state = sim$cell_state,
                                      time = sim$cell_time,
                                      row.names = paste0("c", 1:120)))
  emb <- as.matrix(sim$cell_embedding)
  rownames(emb) <- paste0("c", 1:120)
  colnames(emb) <- c("PC_1", "PC_2")
  dr <- Seurat::CreateDimReducObject(
    embeddings = emb, key = "PC_", assay = "RNA")
  obj[["pca"]] <- dr
  res <- runIOT(obj, group.by = "state", split.by = "time",
                from = "t0", to = "t1", reduction = "pca", root = "S1",
                n_dim = 2)
  expect_equal(dim(res$Q), c(4, 4))
  expect_equal(unname(rowSums(res$Q)), rep(1, 4), tolerance = 1e-8)
})

test_that("runIOT.default raises informative error", {
  expect_error(runIOT(1:3), "embedding matrix")
})
