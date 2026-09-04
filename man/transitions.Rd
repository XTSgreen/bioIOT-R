\name{transition_matrix}
\alias{transition_matrix}
\alias{pseudotime_from_transition}
\alias{build_state_features}
\title{Transition matrix, pseudotime and state features}
\usage{
transition_matrix(fit, which = 1, phi = NULL, a = NULL, b = NULL)
pseudotime_from_transition(Q, root)
build_state_features(u)
}
\arguments{
\item{fit}{A \code{bioIOT_fit} object.}
\item{which}{Scenario index (default 1).}
\item{phi, a, b}{Optional user-supplied scenario overriding the stored
one (\code{phi} must be standardized the same way as in the fit).}
\item{Q}{(K, K) row-stochastic transition matrix.}
\item{root}{Root state: integer index or a name matching rownames of Q.}
\item{u}{(K, D) matrix of state feature vectors (e.g. cluster centroids
in a reduced space).}
}
\value{
\code{transition_matrix}: (K, K) row-stochastic transition matrix.
\code{pseudotime_from_transition}: named numeric vector of expected
random-walk hitting times (root = 0). \code{build_state_features}:
(K, K, D + 1) feature array with dimnames
\code{list(NULL, NULL, c(colnames(u), "sim"))}.
}
\description{
\code{transition_matrix}: re-solves the plan at the fitted \eqn{\theta}
and returns \eqn{Q = P/a}. \code{pseudotime_from_transition}: expected
number of random-walk steps to reach the root state (root absorbing),
turning an IOT transition matrix into per-state pseudotime.
\code{build_state_features}: builds the paper's feature library:
per-dimension pure-column target-state features plus the state-similarity
interaction block \eqn{\phi_{ij,F+1} = u_i \cdot u_j}.
}
\examples{
set.seed(1)
sim <- simulate_iot_states(K = 5, seed = 1)
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
Q <- transition_matrix(fit)
rowSums(Q)
pseudotime_from_transition(Q, root = 1)
u <- matrix(rnorm(6 * 2), 6, 2, dimnames = list(NULL, c("PC1", "PC2")))
dim(build_state_features(u))
}
