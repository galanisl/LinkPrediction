#' Ties between members of a karate club
#'
#' Ties between the members of a karate club as recorded by Wayne Zachary in 
#' 1977.
#'
#' @docType data
#'
#' @format An igraph object representing the network
#'
#' @keywords datasets igraph zachary karate
#'
#' @references Zachary, Wayne (1977) An information flow model for conflict and 
#' fission in small groups. \emph{J. of Anthropological Research} 33:452-473
#'
#' @source \href{http://konect.uni-koblenz.de/networks/ucidata-zachary}{Koblenz 
#' Network Collection.}
#'
#' @examples
#' # Show number of nodes and edges
#' igraph::vcount(karate_club)
#' igraph::ecount(karate_club)
#' 
#' @importFrom igraph vcount ecount
#' 
"karate_club"

