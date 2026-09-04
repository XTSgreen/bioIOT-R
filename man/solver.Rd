\name{soft_sinkhorn}
\alias{soft_sinkhorn}
\alias{row_conditional}
\title{Semi-relaxed OT solver and row-conditional transitions}
\usage{
soft_sinkhorn(C, a, b, mu = 0.5, eps = 1, iters = 1000L, damp = 0.5,
  tol = 1e-12)
row_conditional(P, a, eps_t = 1e-300)
}
\arguments{
\item{C}{(K, K) cost matrix.}
\item{a, b}{(K,) source / target masses; any positive scale is accepted
(auto-normalized). Zero-mass states are allowed in \code{a}.}
\item{mu}{Soft-marginal strength; large \code{mu} approaches hard-marginal
OT, \code{mu = 0} recovers a plain row-softmax.}
\item{eps}{Entropic regularization (> 0).}
\item{iters}{Fixed-point iteration cap.}
\item{damp}{Damping factor of the base iteration.}
\item{tol}{Fixed-point residual tolerance.}
\item{P}{(K, K) transport plan.}
\item{eps_t}{Numerical floor.}
}
\value{
\code{soft_sinkhorn}: (K, K) transport plan with row sums equal to the
normalized \code{a}. \code{row_conditional}: (K, K) matrix; rows with
\eqn{a_i = 0} are all-zero.
}
\description{
\code{soft_sinkhorn} solves the row-marginal-hard / column-marginal-KL-soft
problem
\eqn{\min_P <C,P> - \eps H(P) + \mu KL(col(P) || b)} subject to
\eqn{P 1 = a}, via Anderson-accelerated damped fixed-point iteration.
\code{row_conditional} returns \eqn{Q(i,\cdot) = P(i,\cdot) / a_i}.
}
\examples{
K <- 5
set.seed(1)
P <- soft_sinkhorn(matrix(rnorm(25), 5), rep(1, 5), rep(1, 5))
rowSums(P)
row_conditional(P, rep(1, 5) / 5)
}
