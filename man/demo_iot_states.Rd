\name{demo_iot_states}
\alias{demo_iot_states}
\docType{data}
\title{Bundled bioIOT demo dataset}
\format{
A list as returned by \code{simulate_iot_states(seed = 1)} with default
arguments: \code{u}, \code{phi}, \code{a}, \code{b}, \code{P_true},
\code{T_true}, \code{theta_true}, \code{embedding}, \code{K},
\code{F_emb}, \code{cell_embedding}, \code{cell_state},
\code{cell_time}.
}
\description{
Pre-computed synthetic demo data for quick experimentation with
\code{\link{fit_iot}}, \code{\link{transition_matrix}} and
\code{\link{runIOT}}. Regenerate with
\code{bioIOT::simulate_iot_states(seed = 1)}.
}
\examples{
data(demo_iot_states)
names(demo_iot_states)
demo_iot_states$K
}
\keyword{datasets}
