#' Collaborations between jazz musicians
#'
#' Network of collaborations between jazz musicians. Each node is musician and 
#' they are connected if they have played together.
#'
#' @docType data
#'
#' @format An igraph object representing the network
#'
#' @keywords datasets igraph jazz
#'
#' @references Gleiser P. M. and Danon L. (2003) Community structure in jazz.
#' \emph{Advances in Complex Systems} 6(4):565-573
#'
#' @source \href{http://konect.uni-koblenz.de/networks/arenas-jazz}{Koblenz 
#' Network Collection.}
#'
#' @examples
#' # Show number of nodes and edges
#' igraph::vcount(jazz_collab)
#' igraph::ecount(jazz_collab)
#' 
"jazz_collab"

