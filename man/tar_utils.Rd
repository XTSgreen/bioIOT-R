\name{tar_utils}
\alias{has_cel_file}
\alias{find_platform_file}
\title{Inspect a GEO RAW tar without extracting}
\usage{
has_cel_file(tar_path)
find_platform_file(tar_path)
}
\arguments{
\item{tar_path}{Path to a \code{*_RAW.tar} archive.}
}
\value{
\code{has_cel_file}: TRUE if the tar contains \code{*.CEL} (optionally
\code{.gz}) members. \code{find_platform_file}: the first \code{*.bgx.gz}
platform filename, or \code{NA_character_}.
}
\description{
Check whether a GEO \code{*_RAW.tar} contains CEL data files or only
platform (bgx.gz) files, without extracting the archive.
}
\examples{
\dontrun{
has_cel_file("GSE62254_RAW.tar")
find_platform_file("GSE26901_RAW.tar")
}
}
