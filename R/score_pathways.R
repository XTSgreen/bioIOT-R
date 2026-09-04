# bioIOT pathway scoring
#
# 每基因先跨样本 z-score，再对每个通路的 marker 基因取均值。
# 与 13_bulk_cohort_gex.R 的 score_pathways 完全同源。

#' Score bulk-cohort samples on the 8 IOT pathways
#'
#' Per-gene z-scores across samples, then the mean of the pathway marker
#' genes per sample. Exactly the scoring used to build the paper's bulk
#' pathway scores from gene x sample expression matrices.
#'
#' @param expr_gene Gene x sample numeric matrix (rownames = gene symbols).
#' @param markers Named list of character vectors; default
#'   \code{\link{pathway_markers}} (8 pathways).
#'
#' @return Sample x pathway numeric matrix (rownames = sample ids,
#'   colnames = pathway names). Pathways without any marker present in
#'   \code{expr_gene} yield \code{NA} columns.
#'
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(200 * 6), nrow = 200)
#' rownames(expr) <- paste0("G", 1:200)
#' colnames(expr) <- paste0("S", 1:6)
#' rownames(expr)[1] <- "VIM"; rownames(expr)[2] <- "CDH2"
#' pw <- score_pathways(expr)
#' @export
score_pathways <- function(expr_gene, markers = pathway_markers) {
  if (!is.matrix(expr_gene)) {
    stop("`expr_gene` must be a gene x sample matrix")
  }
  if (nrow(expr_gene) == 0 || ncol(expr_gene) == 0) {
    stop("`expr_gene` must be non-empty")
  }
  z <- t(scale(t(expr_gene)))
  pw_scores <- sapply(names(markers), function(pw) {
    g <- intersect(markers[[pw]], rownames(expr_gene))
    if (length(g) == 0) return(rep(NA_real_, ncol(expr_gene)))
    colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
  rownames(pw_scores) <- colnames(expr_gene)
  pw_scores
}
