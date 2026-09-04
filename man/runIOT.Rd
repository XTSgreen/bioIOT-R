\name{runIOT}
\alias{runIOT}
\alias{runIOT.default}
\alias{runIOT.matrix}
\alias{runIOT.SingleCellExperiment}
\alias{runIOT.Seurat}
\title{Run IOT on a single-cell object}
\usage{
runIOT(object, ...)

\method{runIOT}{matrix}(object, state, from, to, root = NULL,
  n_dim = NULL, T_obs = NULL, mu = 0.5, eps = 1, lam = 0.05,
  epochs = 300, lr = 0.01, n_restart = 4, seed = 1,
  two_stage = TRUE, ...)

\method{runIOT}{SingleCellExperiment}(object, state_col, time_col,
  from, to, dimred = "PCA", root = NULL, ...)

\method{runIOT}{Seurat}(object, group.by = NULL, split.by = NULL,
  from = NULL, to = NULL, reduction = "pca", root = NULL, ...)
}
\arguments{
\item{object}{Cell embedding matrix (cells x dims) or a
\code{SingleCellExperiment} / \code{Seurat} object.}
\item{state}{Factor/character of state (cluster) labels, one per cell.}
\item{from, to}{Logical vectors or indices marking source / target
timepoint cells.}
\item{root}{Optional root state (name or index) for pseudotime.}
\item{n_dim}{Number of embedding dims used as state features.}
\item{T_obs}{Optional (K, K) observed row-conditional transition matrix
(e.g. from lineage/clone data). When given, feature weights are fitted;
otherwise the plan is solved with uniform feature weights.}
\item{mu, eps, lam, epochs, lr, n_restart, seed, two_stage}{Passed to
\code{\link{fit_iot}} (only used when \code{T_obs} is given).}
\item{state_col}{State label column in colData (or a vector).}
\item{time_col}{Timepoint column in colData (or a vector).}
\item{dimred}{Reduced dimension name (default "PCA").}
\item{group.by}{Seurat ident / metadata column for states.}
\item{split.by}{Metadata column holding the timepoint values
\code{from}/\code{to}.}
\item{reduction}{Reduction name (default "pca").}
\item{...}{Additional arguments.}
}
\value{
List with \code{fit} (\code{bioIOT_fit} or NULL), \code{Q} (transition
matrix), \code{u} (state centroids), \code{a}, \code{b}, \code{phi},
\code{theta} and \code{pseudotime} (if \code{root} given).
}
\description{
Aggregates cells into states (clusters), builds state-transition features
from the cell embedding, optionally fits IOT weights against observed
transitions \code{T_obs}, and returns the state transition matrix with
optional random-walk pseudotime.
}
\examples{
set.seed(1)
sim <- simulate_iot_states(K = 5, n_cells = 40, seed = 1)
res <- runIOT(sim$cell_embedding, sim$cell_state,
              from = sim$cell_time == "t0", to = sim$cell_time == "t1",
              root = "S1", n_restart = 1, epochs = 100)
res$Q[1:3, 1:3]
res$pseudotime
}
