\name{pathway_markers}
\alias{pathway_markers}
\title{IOT pathway marker genes (8 pathways)}
\usage{pathway_markers}
\value{
A named list of 8 character vectors of HGNC gene symbols.
}
\description{
The pathway marker gene library used by the paper's 8-dimensional IOT
cost features: EMT, Angiogenesis, Hypoxia, Stemness, Immune_Cytotoxic,
CellCycle, TGFb and Chemokine. Genes are HGNC symbols.
}
\examples{
names(pathway_markers)
length(pathway_markers$EMT)
}
\keyword{datasets}
