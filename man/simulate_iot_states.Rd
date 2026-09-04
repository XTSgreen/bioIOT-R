\name{simulate_iot_states}
\alias{simulate_iot_states}
\title{Simulate single-cell state-transition data for bioIOT}
\usage{
simulate_iot_states(K = 6, F_emb = 2, n_cells = 50, seed = 1,
  theta_true = c(0.9, -0.7, 1.1), mu = 0.5, eps = 1)
}
\arguments{
\item{K}{Number of states.}
\item{F_emb}{Number of state feature dimensions (>= 2 for plotting).}
\item{n_cells}{Cells per state per timepoint.}
\item{seed}{Random seed.}
\item{theta_true}{True feature weights (length F_emb + 1).}
\item{mu}{Solver hyperparameter used to build the true plan.}
\item{eps}{Solver hyperparameter used to build the true plan.}
}
\value{
List with \code{u} (K x F_emb state centroids), \code{phi}
((K, K, F_emb + 1) features), \code{a}, \code{b} (normalized state
masses), \code{P_true}, \code{T_true} (row-conditional),
\code{theta_true}, \code{embedding} (K x 2 state coordinates),
\code{K}, \code{F_emb}, \code{cell_embedding} (cells x F_emb),
\code{cell_state} and \code{cell_time} (cell-level factors).
}
\description{
Generates a reproducible synthetic single-cell-like dataset: states with
feature centroids, source/target timepoint masses, a true IOT plan and
transition matrix, plus cell-level metadata for \code{\link{runIOT}}.
The bundled \code{\link{demo_iot_states}} dataset is produced by this
function with its default seed.
}
\examples{
sim <- simulate_iot_states(K = 5, seed = 1)
names(sim)
round(sim$T_true, 2)
}
