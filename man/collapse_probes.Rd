\name{collapse_probes}
\alias{collapse_probes}
\title{Collapse a probe-level expression matrix to gene level}
\usage{
collapse_probes(expr, probe_id = rownames(expr), gene_id)
}
\arguments{
\item{expr}{Probe x sample numeric matrix.}
\item{probe_id}{Character vector of probe ids, length \code{nrow(expr)}
(defaults to \code{rownames(expr)}).}
\item{gene_id}{Character vector of gene symbols matching \code{probe_id};
NA / empty entries are dropped.}
}
\value{
Gene x sample numeric matrix (rownames = gene symbols).
}
\description{
For each gene symbol keep the probe with the highest cross-sample mean
expression. This is the paper's probe-annotation folding rule.
}
\examples{
expr <- matrix(rnorm(40), nrow = 4)
colnames(expr) <- paste0("S", 1:10)
collapse_probes(expr, c("P1", "P2", "P3", "P4"),
                c("VIM", "VIM", "CDH2", "MYC"))
}
