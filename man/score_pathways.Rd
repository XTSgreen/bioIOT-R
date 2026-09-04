\name{score_pathways}
\alias{score_pathways}
\title{Score bulk-cohort samples on the 8 IOT pathways}
\usage{
score_pathways(expr_gene, markers = pathway_markers)
}
\arguments{
\item{expr_gene}{Gene x sample numeric matrix (rownames = gene symbols).}
\item{markers}{Named list of character vectors; default
\code{\link{pathway_markers}} (8 pathways).}
}
\value{
Sample x pathway numeric matrix (rownames = sample ids, colnames =
pathway names). Pathways without any marker present in \code{expr_gene}
yield \code{NA} columns.
}
\description{
Per-gene z-scores across samples, then the mean of the pathway marker
genes per sample. Exactly the scoring used to build the paper's bulk
pathway scores from gene x sample expression matrices.
}
\examples{
set.seed(1)
expr <- matrix(rnorm(200 * 6), nrow = 200)
rownames(expr) <- paste0("G", 1:200)
colnames(expr) <- paste0("S", 1:6)
rownames(expr)[1] <- "VIM"; rownames(expr)[2] <- "CDH2"
pw <- score_pathways(expr)
dim(pw)
}
