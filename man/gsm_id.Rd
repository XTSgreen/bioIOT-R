\name{gsm_id}
\alias{gsm_id}
\title{Extract GSM ids from CEL filenames}
\usage{
gsm_id(filenames)
}
\arguments{
\item{filenames}{Character vector of filenames such as
\code{"GSM1523727_INT_A.CEL.gz"}.}
}
\value{
Character vector of GSM ids; \code{NA} where no match.
}
\description{
Extract GSM sample ids from CEL filenames, vectorized.
}
\examples{
gsm_id(c("GSM1523727_INT_A.CEL.gz", "sample_no_id.CEL", NA))
}
