\name{make_cost}
\alias{make_cost}
\alias{row_ce_loss}
\alias{zscore_phi}
\title{Cost construction, likelihood and feature standardization}
\usage{
make_cost(phi, theta)
row_ce_loss(T, P, a, eps_t = 1e-12)
zscore_phi(phi)
}
\arguments{
\item{phi}{(K, K, F) feature array; a (K, K) matrix is promoted to a
single feature.}
\item{theta}{Numeric length F (a scalar is recycled for 2-D \code{phi}).}
\item{T}{(K, K) observed row-conditional transition matrix.}
\item{P}{(K, K) model plan.}
\item{a}{(K,) source masses.}
\item{eps_t}{Log floor.}
}
\value{
\code{make_cost}: (K, K) cost matrix. \code{row_ce_loss}: numeric loss
(lower is better). \code{zscore_phi}: list with \code{phi_z}
((K, K, F) standardized array) and \code{meta} (2 x F matrix with rows
"mean" and "sd").
}
\description{
\code{make_cost}: linear feature cost
\eqn{C_{ij} = -\sum_k \phi_{ijk} \theta_k}.
\code{row_ce_loss}: cross-entropy between observed row-conditional
transitions and the model \eqn{Q = P/a}; zero-mass source states are
masked out. \code{zscore_phi}: z-scores each feature across all (i, j)
entries and returns the transform for re-application.
}
\examples{
phi <- array(rnorm(3 * 3 * 2), c(3, 3, 2))
C <- make_cost(phi, c(1, -1))
out <- zscore_phi(phi)
row_ce_loss(matrix(0.5, 3, 3), matrix(0.25, 3, 3), rep(1, 3) / 3)
}
