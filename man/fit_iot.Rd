\name{fit_iot}
\alias{fit_iot}
\alias{print.bioIOT_fit}
\alias{summary.bioIOT_fit}
\title{Fit inverse optimal transport feature weights}
\usage{
fit_iot(phi, a, b, T, mu = 0.5, lam = 0.05, eps = 1, epochs = 300,
  lr = 0.01, n_restart = 4, seed = 1, two_stage = TRUE,
  iters = 1000L, damp = 0.5, standardize = TRUE, verbose = FALSE)
\method{print}{bioIOT_fit}(x, ...)
\method{summary}{bioIOT_fit}(object, ...)
}
\arguments{
\item{phi}{(K, K, F) feature array, or a list of arrays (scenarios).}
\item{a, b}{(K,) masses, or lists (auto-normalized).}
\item{T}{(K, K) observed row-conditional transitions, or a list.}
\item{mu}{Soft-marginal strength (default 0.5, paper working point).}
\item{lam}{l1 strength in stage 1 (default 0.05).}
\item{eps}{Entropic regularization (default 1).}
\item{epochs}{Adam steps per stage (default 300).}
\item{lr}{Adam learning rate (default 0.01).}
\item{n_restart}{Random restarts (default 4).}
\item{seed}{Base seed for restart inits.}
\item{two_stage}{Debias refit on the selected support (default TRUE).}
\item{iters, damp}{Forward solver controls.}
\item{standardize}{Z-score each scenario's features before fitting
(default TRUE, paper pipeline).}
\item{verbose}{Print progress.}
\item{x}{A \code{bioIOT_fit} object.}
\item{object}{A \code{bioIOT_fit} object.}
\item{...}{Unused.}
}
\value{
Object of class \code{bioIOT_fit} with components \code{theta}
(debiased weights), \code{theta1} (stage-1 weights), \code{support}
(logical mask), \code{loss} (final row-CE), \code{restart_losses},
\code{scenarios} (stored standardized scenarios), \code{meta}
(standardization transforms), \code{hyper} and dims \code{K}/\code{F}.
}
\description{
Fits feature weights \eqn{\theta} so that the soft-marginal OT plan
induced by \eqn{C = -einsum(\phi, \theta)} reproduces the observed
row-conditional transitions \eqn{T}. Two-stage fitting (l1 selection then
Buehlmann-style debias refit) with multi-restart; gradients are exact
implicit differentiations of the fixed point (implicit function theorem),
which stays numerically stable where unrolled backpropagation diverges.
}
\examples{
set.seed(1)
sim <- simulate_iot_states(K = 5, seed = 1)
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
fit$theta
summary(fit)
}
