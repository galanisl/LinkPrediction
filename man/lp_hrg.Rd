% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/link_predictors.R
\name{lp_hrg}
\alias{lp_hrg}
\title{Link prediction with the HRG model}
\usage{
lp_hrg(g, samples = 1000)
}
\arguments{
\item{g}{igraph; The network of interest.}

\item{samples}{integer; Number of HRG models to consider.}
}
\value{
Tibble with the following columns:
\item{nodeA}{The ID of a network node.}
\item{nodeB}{The ID of a network node.}
\item{scr}{The likelihood score of interaction for the node pair.}
}
\description{
Given a network of interest, it computes the likelihood score of interaction,
for all disconnected node pairs, based on a Hierarchical Random Graph model
(HRG).
}
\examples{
# Apply the HRG link predictor to the Zachary Karate Club network
hrg <- lp_hrg(g = karate_club, samples = 100)

}
\references{
Clauset, A. et al. (2008) Hierarchical structure and the 
prediction of missing links in networks. \emph{Nature} 453(7191):98-101
}
\author{
Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
}
