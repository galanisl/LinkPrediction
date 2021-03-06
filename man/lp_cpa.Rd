% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/link_predictors.R
\name{lp_cpa}
\alias{lp_cpa}
\title{Link prediction with CPA}
\usage{
lp_cpa(g)
}
\arguments{
\item{g}{igraph; The network of interest.}
}
\value{
Tibble with the following columns:
\item{nodeA}{The ID of a network node.}
\item{nodeB}{The ID of a network node.}
\item{scr}{The likelihood score of interaction for the node pair.}
}
\description{
Given a network of interest, it computes the likelihood score of interaction,
for all disconnected node pairs, with the CAR-based Preferential Attachment
(CPA)
}
\examples{
# Apply the CPA link predictor to the Zachary Karate Club network
cpa <- lp_cpa(g = karate_club)

}
\references{
Cannistraci, C. V. et al. (2013) From link-prediction in brain 
connectomes and protein interactomes to the local-community-paradigm in 
complex networks. \emph{Sci. Rep.} 3:1613
}
\author{
Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
}
