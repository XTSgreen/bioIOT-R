# bioIOT probe collapse
#
# 同一 gene symbol 对应多个探针时，保留跨样本平均表达最高的探针
# （与 13_bulk_cohort_gex.R 的注释折叠逻辑同源）。

#' Collapse a probe-level expression matrix to gene level
#'
#' For each gene symbol keep the probe with the highest cross-sample mean
#' expression. This is the paper's probe-annotation folding rule.
#'
#' @param expr Probe x sample numeric matrix.
#' @param probe_id Character vector of probe ids, length nrow(expr)
#'   (defaults to rownames(expr)).
#' @param gene_id Character vector of gene symbols matching probe_id
#'   (NA / "" entries are dropped).
#'
#' @return Gene x sample numeric matrix (rownames = gene symbols).
#'
#' @examples
#' expr <- matrix(rnorm(40), nrow = 4)
#' colnames(expr) <- paste0("S", 1:10)
#' probe_id <- c("P1", "P2", "P3", "P4")
#' gene_id <- c("VIM", "VIM", "CDH2", "MYC")
#' collapse_probes(expr, probe_id, gene_id)
#' @export
collapse_probes <- function(expr, probe_id = rownames(expr), gene_id) {
  if (!is.matrix(expr)) stop("`expr` must be a probe x sample matrix")
  if (length(probe_id) != nrow(expr)) {
    stop("`probe_id` length must equal nrow(expr)")
  }
  if (length(gene_id) != length(probe_id)) {
    stop("`gene_id` length must equal `probe_id` length")
  }
  probe_id <- as.character(probe_id)
  gene_id <- as.character(gene_id)
  keep <- !is.na(gene_id) & nzchar(gene_id)
  if (!any(keep)) stop("no valid gene symbols")
  expr <- expr[keep, , drop = FALSE]
  probe_id <- probe_id[keep]
  gene_id <- gene_id[keep]

  probe_mean <- rowMeans(expr, na.rm = TRUE)
  ord <- order(gene_id, -probe_mean)
  expr <- expr[ord, , drop = FALSE]
  gene_id <- gene_id[ord]
  first <- !duplicated(gene_id)
  out <- expr[first, , drop = FALSE]
  rownames(out) <- gene_id[first]
  out
}
