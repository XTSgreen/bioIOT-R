\name{plot_transition_heatmap}
\alias{plot_transition_heatmap}
\alias{plot_transition_flow}
\alias{plot_theta}
\alias{plot_pathway_trend}
\title{bioIOT visualisation}
\usage{
plot_transition_heatmap(Q, labels = NULL)
plot_transition_flow(Q, embedding, labels = NULL, threshold = 0.05)
plot_theta(fit, labels = NULL)
plot_pathway_trend(pw_scores, time, se = TRUE)
}
\arguments{
\item{Q}{(K, K) row-stochastic transition matrix.}
\item{labels}{Optional state / feature labels.}
\item{embedding}{(K x >= 2) state coordinates (e.g. cluster centroids).}
\item{threshold}{Minimum transition mass to draw an arrow.}
\item{fit}{A \code{bioIOT_fit} object (or named numeric weights with an
optional \code{support} attribute).}
\item{pw_scores}{Sample x pathway score matrix.}
\item{time}{Numeric time / pseudotime per sample (row).}
\item{se}{Show the smooth confidence band (default TRUE).}
}
\value{
A ggplot object.
}
\description{
\code{plot_transition_heatmap}: transition-matrix heatmap with annotated
weights. \code{plot_transition_flow}: transition arrows on a 2-D state
embedding, width proportional to transition mass (CellRank-style).
\code{plot_theta}: barplot of fitted feature weights, selected support
highlighted. \code{plot_pathway_trend}: pathway score trajectories over
time (loess).
}
\examples{
Q <- matrix(c(0.7, 0.3, 0, 0.2, 0.6, 0.2, 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
plot_transition_heatmap(Q)
emb <- matrix(c(0, 0, 1, 1, 2, 0), 3, 2, byrow = TRUE)
plot_transition_flow(Q, emb)

set.seed(1)
sim <- simulate_iot_states(K = 5, seed = 1)
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
plot_theta(fit)

expr <- matrix(rnorm(200 * 12), nrow = 200)
rownames(expr) <- paste0("G", 1:200)
colnames(expr) <- paste0("S", 1:12)
rownames(expr)[1] <- "VIM"; rownames(expr)[2] <- "CDH2"; rownames(expr)[3] <- "MKI67"
plot_pathway_trend(score_pathways(expr), time = sort(runif(12)))
}
