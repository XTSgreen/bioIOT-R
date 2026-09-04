test_that("pathway_markers: 8 pathways with valid symbols", {
  expect_length(pathway_markers, 8)
  expect_named(pathway_markers,
               c("EMT", "Angiogenesis", "Hypoxia", "Stemness",
                 "Immune_Cytotoxic", "CellCycle", "TGFb", "Chemokine"))
  ok <- vapply(pathway_markers, function(g)
    all(grepl("^[A-Z][A-Z0-9-]*$", g)), logical(1))
  expect_true(all(ok))
})

test_that("score_pathways: z-score then marker mean", {
  set.seed(1)
  n_gene <- 200
  expr <- matrix(rnorm(n_gene * 6), nrow = n_gene)
  rownames(expr) <- paste0("G", seq_len(n_gene))
  colnames(expr) <- paste0("S", 1:6)
  rownames(expr)[1] <- "VIM"
  rownames(expr)[2] <- "CDH2"
  pw <- score_pathways(expr)
  expect_equal(dim(pw), c(6, 8))
  expect_equal(rownames(pw), colnames(expr))
  # custom single-pathway library: mean of z-scored markers
  m <- list(Two = c("VIM", "CDH2"))
  z <- t(scale(t(expr)))
  ref <- colMeans(z[c("VIM", "CDH2"), , drop = FALSE])
  got <- score_pathways(expr, m)[, "Two"]
  expect_equal(unname(got), unname(ref))
  # absent pathway -> NA column
  m2 <- list(Nope = c("ZZZ1"))
  pw2 <- score_pathways(expr, m2)
  expect_true(all(is.na(pw2[, "Nope"])))
})

test_that("collapse_probes: keeps highest-mean probe per gene", {
  expr <- matrix(c(rep(1, 10), rep(5, 10), rep(2, 10), rep(3, 10)),
                 nrow = 4, byrow = TRUE)
  colnames(expr) <- paste0("S", 1:10)
  out <- collapse_probes(expr, c("P1", "P2", "P3", "P4"),
                         c("VIM", "VIM", "CDH2", "MYC"))
  expect_equal(nrow(out), 3)
  expect_equal(rownames(out), c("CDH2", "MYC", "VIM"))   # alphabetical
  expect_equal(out["VIM", 1], 5)   # P2 kept over P1
  expect_equal(out["CDH2", 1], 2)
})

test_that("gsm_id: extracts ids, NA otherwise", {
  got <- gsm_id(c("GSM1523727_INT_A.CEL.gz", "sample_no_id.CEL", NA))
  expect_equal(got[1], "GSM1523727")
  expect_true(is.na(got[2]))
  expect_true(is.na(got[3]))
})

test_that("input validation raises informative errors", {
  expect_error(score_pathways(matrix(0, 0, 3)), "non-empty")
  expect_error(collapse_probes(matrix(1:4, 2), gene_id = c("A")),
               "length must equal")
  expect_true(all(is.na(gsm_id("no_id.CEL"))))
})
