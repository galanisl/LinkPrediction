
globalVariables(c("nodeA", "nodeB", "scr", "spm"))

#' Disconnected node pairs in a network
#' 
#' Given a network of interest with N nodes and L links, produces a matrix with
#' N(N-1)/2 - L rows and 2 columns containing the list of disconnected node
#' pairs in the network.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return A matrix with the disconnected node pairs in the network.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Determine which node pairs are disconnected in the Zacharay Karate Club 
#' # network
#' non_edges <- get_non_edges(g = karate_club)
#' 
#' @export
#' @importFrom igraph as_adjacency_matrix
#' @import Matrix
#'
get_non_edges <- function(g){
  m <- as_adjacency_matrix(g, names = FALSE)
  m[lower.tri(m, diag = TRUE)] <- 1
  return(which(m == 0, arr.ind = TRUE))
}


# Neighbourhood-based link predictors -------------------------------------

#' Link prediction with Common Neighbours
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Common Neighbours link predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Newman, M. E. J. (2001) Clustering and preferential attachment in
#' growing networks. \emph{Phys. Rev. E.} 64:025102
#'
#' @examples
#' # Apply the Common Neighbours link predictor to the Zachary Karate Club 
#' # network
#' cn <- lp_cn(g = karate_club)
#' 
#' @export
#' @importFrom igraph as_adjacency_matrix
#' @import Matrix
#' @importFrom dplyr tibble %>% arrange desc
#'
lp_cn <- function(g){
  non_edges <- get_non_edges(g)
  m <- as_adjacency_matrix(g, names = FALSE)
  cn <- m %*% m
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = cn[non_edges]) %>% 
    arrange(desc(scr))
  return(prediction)
}

#' Link prediction with Preferential Attachment
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Preferential Attachment link 
#' predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Newman, M. E. J. (2001) Clustering and preferential attachment in
#' growing networks. \emph{Phys. Rev. E.} 64:025102
#'
#' @examples
#' # Apply the Preferential Attachment link predictor to the Zachary Karate Club 
#' # network
#' pa <- lp_pa(g = karate_club)
#' 
#' @export
#' @importFrom igraph degree
#' @importFrom dplyr tibble %>% arrange desc
#'
lp_pa <- function(g){
  non_edges <- get_non_edges(g)
  pa <- degree(g) %*% t(degree(g))
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = pa[non_edges]) %>% 
    arrange(desc(scr))
  return(prediction)
}

#' Link prediction with methods available in igraph
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the link predictors provided by igraph.
#' This function is used by \code{lp_jc}, \code{lp_dice} and \code{lp_aa}.
#' 
#' @param g igraph; The network of interest.
#' @param method character; One of "jaccard", "dice" or "invlogweighted".
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @importFrom igraph similarity
#' @import Matrix
#' @importFrom dplyr tibble %>% arrange desc
#'
lp_igraph <- function(g, method){
  non_edges <- get_non_edges(g)
  sim <- similarity(g, method = method)
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = sim[non_edges]) %>% 
    arrange(desc(scr))
  return(prediction)
}

#' Link prediction with Jaccard's index
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Jaccard link predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Jaccard, P. (1912) The distribution of the flora in the alpine 
#' zone. \emph{New phytologist} 11(2):37-50
#'
#' @examples
#' # Apply the Jaccard link predictor to the Zachary Karate Club network
#' jc <- lp_jc(g = karate_club)
#' 
#' @export
#'
lp_jc <- function(g){
  prediction <- lp_igraph(g, method = "jaccard")
  return(prediction)
}

#' Link prediction with Dice's index
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Dice link predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Dice, L. R. (1945) Measures of the amount of ecologic 
#' association between species. \emph{Ecology} 26(3):297-302
#'
#' @examples
#' # Apply the Dice link predictor to the Zachary Karate Club network
#' dice <- lp_dice(g = karate_club)
#' 
#' @export
#'
lp_dice <- function(g){
  prediction <- lp_igraph(g, method = "dice")
  return(prediction)
}

#' Link prediction with Adamic and Adar's index
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Adamic and Adar link predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Adamic, L. and Adar, E. (2003) Friends and neighbors on the Web. 
#' \emph{Social Networks} 25(3):211-230
#'
#' @examples
#' # Apply the Adamic and Adar link predictor to the Zachary Karate Club network
#' aa <- lp_aa(g = karate_club)
#' 
#' @export
#'
lp_aa <- function(g){
  prediction <- lp_igraph(g, method = "invlogweighted")
  return(prediction)
}

#' Link prediction with Resource Allocation
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Resource Allocation link predictor.
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Zhou, T. et al. (2009) Predicting missing links via local 
#' information. \emph{The European Physical Journal B} 71(4):623-630
#'
#' @examples
#' # Apply the Resource Allocation link predictor to the Zachary Karate Club 
#' # network
#' ra <- lp_ra(g = karate_club)
#' 
#' @export
#' @importFrom igraph adjacent_vertices V
#' @importFrom dplyr tibble %>% arrange desc
#' @importFrom purrr map2_dbl
#'
lp_ra <- function(g){
  non_edges <- get_non_edges(g)
  # Compute the neighbourhoods of all nodes
  ne <- adjacent_vertices(g, v = V(g))
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = map2_dbl(as.list(nodeA), as.list(nodeB), 
                                         compute_ra, 
                                         ne = ne)) %>% arrange(desc(scr))
  return(prediction)
}

#' Compute the Resource Allocation score of two nodes
#' 
#' Given two indices representing two nodes and a list of neighbourhoods,
#' compute the Resource Allocation index for the node pair. This function is 
#' used by \code{lp_ra}.
#' 
#' @param ne list; The neighbourhood of all nodes from a network of interest.
#' @param i integer; The ID of a network node.
#' @param j integer; The ID of a network node.
#' 
#' @return The Resource Allocation likelihood score of interaction for the node
#' pair i-j.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Zhou, T. et al. (2009) Predicting missing links via local 
#' information. \emph{The European Physical Journal B} 71(4):623-630
#'
#'
compute_ra <- function(ne, i, j){
  # The common neighbourhood of the nodes i and j
  cn <- intersect(ne[[i]], ne[[j]])
  
  # Degrees of the common neighbours
  deg <- lengths(ne[cn])
  
  # Compute the RA index
  ra <- sum(1/deg)
  
  return(ra)
}

# (Dis-)similarity matrix -------------------------------------------------

#' Link prediction with pairwise node (dis-)similarities 
#' 
#' Given a network of interest and a matrix of pairwise node (dis-)similarities
#' produced by some algorithm, it computes the likelihood score of interaction,
#' for all disconnected node pairs, based on such (dis-)similarities.
#' 
#' @param g igraph; The network of interest.
#' @param m matrix; A square matrix of pairwise node (dis-)similarities.
#' @param type character; One of "sim" or "dis", which indicate whether the 
#' matrix stores similarities (high values -> higher likelihood of interaction) 
#' or dissimilarities (low values -> higher likelihood of interaction).
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Compute a node similarity matrix with one of the algorithms included in 
#' # igraph
#' dice <- igraph::similarity(karate_club, method = "dice")
#' 
#' # Now compute likelihood scores for disconnected node pairs with this matrix
#' sim <- lp_matrix(g = karate_club, dice, type = "sim")
#' 
#' @export
#' @import Matrix
#' @importFrom dplyr tibble arrange desc
#'
lp_matrix <- function(g, m, type = "sim"){
  if(!(type %in% c("sim", "dis"))){
    stop("Invalid type provided, valid types are 'sim' and 'dis'")
  }
  non_edges <- get_non_edges(g)
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = m[non_edges])
  
  if(type == "sim"){
    prediction <- arrange(prediction, desc(scr))
  }else if(type == "dis"){
    prediction <- arrange(prediction, scr)
  }
  
  return(prediction)
}

# Embedding-based link predictors -----------------------------------------

#' Link prediction with ISOMAP
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, based on the embedding of the network to 
#' d-dimensional Euclidean space with ISOMAP.
#' 
#' @param g igraph; The network of interest.
#' @param d integer; The dimension of the embedding space.
#' @param use_weights logical; Indicates whether edge weights should be used to
#' compute shortest paths between nodes.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Tenenbaum, J. B. (2000)  A global geometric framework for 
#' nonlinear dimensionality reduction. \emph{Science} 290:2319-2323
#' @references Kuchaiev, O. et al. (2009)  Geometric de-noising of 
#' protein-protein interaction networks. \emph{PLoS Comput. Biol.} 5:e1000454
#'
#' @examples
#' # Apply the ISOMAP link predictor to the Zachary Karate Club network
#' iso <- lp_isomap(g = karate_club, d = 2, use_weights = FALSE)
#' 
#' @export
#' @importFrom igraph distances
#' @import Matrix
#' @importFrom RSpectra svds
#' @importFrom stats dist
#'
lp_isomap <- function(g, d = 2, use_weights = FALSE){
  # Kernel computation
  if(use_weights){
    sp <- distances(g)  
  }else{
    sp <- distances(g, weights = NA)  
  }
  
  # Kernel centring
  N <- nrow(sp)
  J <- diag(N) - (1/N)*matrix(1, N, N) #Form the centring matrix J
  sp <- (-0.5)*(J %*% (sp^2) %*% J)
  
  # Embedding to d-dimensional space
  res <- svds(sp, k = d)
  L <- diag(res$d)
  V <- res$v
  
  node_coords <- t(sqrt(L) %*% t(V))
  
  # Pairwise distances between nodes
  D <- as.matrix(dist(x = node_coords, method = "euclidean"))
  
  # Predict edges based on node distances in the embedding space
  prediction <- lp_matrix(g, m = D, type = "dis")
  
  return(prediction)
}

#' Link prediction with Laplacian Eigenmaps
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, based on the embedding of the network to 
#' d-dimensional Euclidean space with Laplacian Eigenmaps.
#' 
#' @param g igraph; The network of interest.
#' @param d integer; The dimension of the embedding space.
#' @param use_weights logical; Indicates whether edge weights should be used to
#' compute the network's Laplacian matrix.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references  Belkin, M. and Niyogi, P. (2001) Laplacian eigenmaps and 
#' spectral techniques for embedding and clustering. \emph{Adv. Neur. In.} 
#' 14:585-591
#'
#' @examples
#' # Apply the Laplacian Eigenmaps link predictor to the Zachary Karate Club 
#' # network
#' leig <- lp_leig(g = karate_club, d = 2, use_weights = FALSE)
#' 
#' @export
#' @importFrom igraph laplacian_matrix
#' @import Matrix
#' @importFrom RSpectra eigs_sym
#' @importFrom stats dist
#'
lp_leig <- function(g, d = 2, use_weights = FALSE){
  # Laplacian matrix computation
  if(use_weights){
    L <- laplacian_matrix(g)  
  }else{
    L <- laplacian_matrix(g, weights = NA)  
  }
  
  # Embedding to d-dimensional space
  leig <- eigs_sym(L, k = d + 1, which = "LM", sigma = 0)
  
  # Pairwise distances between nodes
  D <- as.matrix(dist(x = leig$vectors[, 1:d], method = "euclidean"))
  
  # Predict edges based on node distances in the embedding space
  prediction <- lp_matrix(g, m = D, type = "dis")
  
  return(prediction)
}

#' Link prediction with MCE or ncMCE
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, based on the embedding of the network to 
#' d-dimensional Euclidean space with centred or non-centred Minimum Curvilinear
#' Embedding (MCE and ncMCE, respectively).
#' 
#' @param g igraph; The network of interest.
#' @param d integer; The dimension of the embedding space.
#' @param centre logical; Indicates whether the Minimum Curvilinear kernel 
#' should be centred or not.
#' @param use_weights logical; Indicates whether edge weights should be used to
#' compute the network's Minimum Curvilinear kernel.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) Minimum curvilinearity to 
#' enhance topological prediction of protein interactions by network embedding.
#' \emph{Bioinformatics} 29:i199-i209 
#'
#' @examples
#' # Apply the ncMCE link predictor to the Zachary Karate Club 
#' # network
#' mce <- lp_mce(g = karate_club, d = 2, centre = FALSE, use_weights = FALSE)
#' 
#' @export
#' @importFrom igraph mst distances
#' @import Matrix
#' @importFrom RSpectra svds
#' @importFrom stats dist
#'
lp_mce <- function(g, d = 2, centre = FALSE, use_weights = FALSE){
  # Kernel computation
  if(use_weights){
    g_mst <- mst(g)  
  }else{
    g_mst <- mst(g, algorithm = "unweighted") 
  }
  
  # Pairwise shortest-paths over the MST
  kernel <- distances(g_mst)
  
  # Kernel centring
  if(centre){
    N <- nrow(kernel)
    J <- diag(N) - (1 / N) * matrix(1, N, N) #Form the centring matrix J
    kernel <- (-0.5)*(J %*% (kernel^2) %*% J)
  }
  
  # Embedding to d-dimensional space
  res <- svds(kernel, k = d)
  L <- diag(res$d)
  V <- res$v
  
  node_coords <- t(sqrt(L) %*% t(V))
  
  # Pairwise distances between nodes
  D <- as.matrix(dist(x = node_coords, method = "euclidean"))
  
  # Predict edges based on node distances in the embedding space
  prediction <- lp_matrix(g, m = D, type = "dis")
  
  return(prediction)
}

# HRG ---------------------------------------------------------------------

#' Link prediction with the HRG model
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, based on a Hierarchical Random Graph model
#' (HRG).
#' 
#' @param g igraph; The network of interest.
#' @param samples integer; Number of HRG models to consider.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Clauset, A. et al. (2008) Hierarchical structure and the 
#' prediction of missing links in networks. \emph{Nature} 453(7191):98-101
#'
#' @examples
#' # Apply the HRG link predictor to the Zachary Karate Club network
#' hrg <- lp_hrg(g = karate_club, samples = 100)
#' 
#' @export
#' @importFrom igraph predict_edges
#' @importFrom dplyr tibble %>% arrange desc
#'
lp_hrg <- function(g, samples = 1000){
   hrg_pred <- predict_edges(g, num.samples = samples)
   prediction <- tibble(nodeA = hrg_pred$edges[, 1], 
                        nodeB = hrg_pred$edges[, 2],
                        scr = hrg_pred$prob) %>% 
     arrange(desc(scr))

   return(prediction)
}

# CAR-based link predictors -----------------------------------------------

#' Link prediction with the CAR index
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Cannistraci-Alanis-Ravasi index 
#' (CAR).
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @examples
#' # Apply the CAR link predictor to the Zachary Karate Club network
#' car <- lp_car(g = karate_club)
#' 
#' @export
#' @importFrom igraph adjacent_vertices V
#' @importFrom dplyr tibble %>% arrange desc
#' @importFrom purrr map2_dbl
#'
lp_car <- function(g){
  non_edges <- get_non_edges(g)
  
  # Compute the neighbourhoods of all nodes
  ne <- adjacent_vertices(g, v = V(g))
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = map2_dbl(as.list(nodeA), as.list(nodeB), 
                                         compute_car, ne = ne)) %>% 
    arrange(desc(scr))
  
  return(prediction)
}

#' Compute the CAR score of two nodes
#' 
#' Given two indices representing two nodes and a list of neighbourhoods,
#' compute the CAR index for the node pair. This function is used by 
#' \code{lp_car}.
#' 
#' @param ne list; The neighbourhood of all nodes from a network of interest.
#' @param i integer; The ID of a network node.
#' @param j integer; The ID of a network node.
#' 
#' @return The CAR likelihood score of interaction for the node pair i-j.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @importFrom purrr reduce
#'
compute_car <- function(ne, i, j){
  # The common neighbourhood of the nodes i and j
  cn <- base::intersect(ne[[i]], ne[[j]])
  
  # Identify any local community (i.e. links between common neighbours)
  lc <- base::intersect(reduce(ne[cn], base::union, .init = NULL), cn)
  
  # Compute the CAR index
  car <- length(cn) * ifelse(length(lc) == 0, 0, length(lc) - 1)
  
  return(car)
}

#' Link prediction with CPA
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the CAR-based Preferential Attachment
#' (CPA)
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @examples
#' # Apply the CPA link predictor to the Zachary Karate Club network
#' cpa <- lp_cpa(g = karate_club)
#' 
#' @export
#' @importFrom igraph adjacent_vertices V
#' @importFrom dplyr tibble %>% arrange desc
#' @importFrom purrr map2_dbl
#'
lp_cpa <- function(g){
  non_edges <- get_non_edges(g)
  
  # Compute the neighbourhoods and degrees of all nodes
  ne <- adjacent_vertices(g, v = V(g))
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = map2_dbl(as.list(nodeA), as.list(nodeB), 
                                         compute_cpa, ne = ne)) %>% 
    arrange(desc(scr))
  
  return(prediction)
}

#' Compute the CPA score of two nodes
#' 
#' Given two indices representing two nodes and a list of neighbourhoods,
#' compute the CPA index for the node pair. This function is used by 
#' \code{lp_cpa}.
#' 
#' @param ne list; The neighbourhood of all nodes from a network of interest.
#' @param i integer; The ID of a network node.
#' @param j integer; The ID of a network node.
#' 
#' @return The CPA likelihood score of interaction for the node pair i-j.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @importFrom purrr reduce
#'
compute_cpa <- function(ne, i, j){
  # The common neighbourhood of the nodes i and j
  cn <- base::intersect(ne[[i]], ne[[j]])
  
  # Identify any local community (i.e. links between common neighbours)
  lc <- base::intersect(reduce(ne[cn], base::union, .init = NULL), cn)
  
  # Compute number of links going outside the local community
  extI <- length(ne[[i]]) - length(cn)
  extJ <- length(ne[[j]]) - length(cn)
  
  # Compute the CAR index
  car <- length(cn) * ifelse(length(lc) == 0, 0, length(lc) - 1)
  
  # Compute the CPA index
  cpa <- (extI * extJ) + (extI * car) + (extJ * car) + car^2
  
  return(cpa)
}

#' Link prediction with CAA
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with CAR-based Adamic and Adar (CAA).
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @examples
#' # Apply the CAA link predictor to the Zachary Karate Club network
#' caa <- lp_caa(g = karate_club)
#' 
#' @export
#' @importFrom igraph adjacent_vertices V
#' @importFrom dplyr tibble %>% arrange desc
#' @importFrom purrr map2_dbl
#'
lp_caa <- function(g){
  non_edges <- get_non_edges(g)
  
  # Compute the neighbourhoods and degrees of all nodes
  ne <- adjacent_vertices(g, v = V(g))
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = map2_dbl(as.list(nodeA), as.list(nodeB), 
                                         compute_caa_ra, ne = ne, 
                                         type = "CAA")) %>% 
    arrange(desc(scr))
  
  return(prediction)
}

#' Link prediction with CRA
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with CAR-based Resource Allocation (CRA).
#' 
#' @param g igraph; The network of interest.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @examples
#' # Apply the CRA link predictor to the Zachary Karate Club network
#' cra <- lp_cra(g = karate_club)
#' 
#' @export
#' @importFrom igraph adjacent_vertices V
#' @importFrom dplyr tibble %>% arrange desc
#' @importFrom purrr map2_dbl
#'
lp_cra <- function(g){
  non_edges <- get_non_edges(g)
  
  # Compute the neighbourhoods and degrees of all nodes
  ne <- adjacent_vertices(g, v = V(g))
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = map2_dbl(as.list(nodeA), as.list(nodeB), 
                                         compute_caa_ra, ne = ne, 
                                         type = "CRA")) %>% 
    arrange(desc(scr))
  
  return(prediction)
}

#' Compute the CAA or CAA score of two nodes
#' 
#' Given two indices representing two nodes and a list of neighbourhoods,
#' compute the CAA or CRA index for the node pair. This function is used by 
#' \code{lp_caa} and \code{lp_cra}.
#' 
#' @param ne list; The neighbourhood of all nodes from a network of interest.
#' @param i integer; The ID of a network node.
#' @param j integer; The ID of a network node.
#' @param type character; One of "CAA" or "CRA".
#' 
#' @return The CAA or CRA likelihood score of interaction for the node pair i-j.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Cannistraci, C. V. et al. (2013) From link-prediction in brain 
#' connectomes and protein interactomes to the local-community-paradigm in 
#' complex networks. \emph{Sci. Rep.} 3:1613
#'
#' @importFrom purrr map
#'
compute_caa_ra <- function(ne, i, j, type = "CAA"){
  # The common neighbourhood of the nodes i and j
  cn <- base::intersect(ne[[i]], ne[[j]])
  
  # Degrees of the common neighbours
  deg <- lengths(ne[cn])
  
  # Degrees of the common neightbours inside the local community
  deg_lc <- lengths(map(ne[cn], base::intersect, cn))
  
  # Compute the CAA or CRA index
  if(type == "CAA"){
    scr <- sum(deg_lc/log2(deg))  
  }else if(type == "CRA"){
    scr <- sum(deg_lc/deg)
  }
  
  return(scr)
}

# Structural Perturbation Method ------------------------------------------

#' Link prediction with the SPM
#' 
#' Given a network of interest, it computes the likelihood score of interaction,
#' for all disconnected node pairs, with the Structural Perturbation Method 
#' (SPM).
#' 
#' @param g igraph; The network of interest.
#' @param p_H numeric; Fraction of network links to remove (perturbation set).
#' @param epochs integer; Number of perturbation sets to consider.
#' @param k integer; If k != NA and k < N (number of network nodes), the method 
#' is approximated by computing the k larget eigenvalues instead of all N.
#' 
#' @return Tibble with the following columns:
#' \item{nodeA}{The ID of a network node.}
#' \item{nodeB}{The ID of a network node.}
#' \item{scr}{The likelihood score of interaction for the node pair.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Lu, L. et al. (2015) Toward link predictability of complex 
#' networks. \emph{PNAS} 112(8):2325-2330
#'
#' @examples
#' # Apply the approximated SPM link predictor to the Zachary Karate Club 
#' # network
#' spm <- lp_spm(g = karate_club, p_H = 0.1, epochs = 10, 
#'               k = round(sqrt(igraph::vcount(karate_club))))
#' 
#' @export
#' @importFrom igraph ecount vcount as_adjacency_matrix delete_edges E E<-
#' @import Matrix
#' @importFrom RSpectra eigs_sym
#' @importFrom dplyr tibble %>% arrange desc
#'
lp_spm <- function(g, p_H = 0.1, epochs = 10, k = NA){
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  L <- ecount(g)
  N <- vcount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  
  # Allocate space for the perturbed matrix
  A_perturbed <- Matrix(0, N, N)
  
  for(i in 1:epochs){
    # Determine the link perturbation set dlt_E at random
    to_remove <- sample(L, links_to_remove)
    g_perturbed <- delete_edges(g, to_remove)
    
    # Compute the adjacency matrix of the graph after removal of dlt_E
    A_R <- as_adjacency_matrix(g_perturbed, type = "both")
    
    # Diagonalise A_R
    if(is.na(k) || k == N){
      XLX <- eigen(A_R, symmetric = T)
    }else{
      XLX <- eigs_sym(A_R, k = k, which = "LA")
    }
    
    # Compute the adjacency matrix of the perturbation set
    dlt_A <- A - A_R
    
    # Compute correction terms for eigenvalues based on dlt_A
    dlt_l <- rowSums((t(XLX$vectors) %*% dlt_A) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # Perturb A_R but keep eigenvectors unchanged
    A_perturbed <- A_perturbed + 
      (XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors))
  }
  # Compute the average of the perturbed matrix
  # If A is highly regular, the random removal dlt_E should not change its 
  # structure too much and A and A_perturbed should be close to each other
  A_perturbed <- A_perturbed / epochs
  
  non_edges <- get_non_edges(g)
  
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = A_perturbed[non_edges]) %>% arrange(desc(scr))
  return(prediction)
}

#' Structural consistency of a network
#' 
#' Given a network of interest, it computes its Structural Consistency sigma_c.
#' High values of sigma_c correlate with higher predictability of missing links
#' in a network.
#' 
#' @param g igraph; The network of interest.
#' @param p_H numeric; Fraction of network links to remove (perturbation set).
#' @param epochs integer; Number of perturbation sets to consider.
#' @param k integer; If k != NA and k < N (number of network nodes), the method 
#' is approximated by computing the k larget eigenvalues instead of all N.
#' 
#' @return A vector with \code{epochs} elements containing the consistency 
#' values for each perturbation set.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Lu, L. et al. (2015) Toward link predictability of complex 
#' networks. \emph{PNAS} 112(8):2325-2330
#'
#' @examples
#' # Compute the structural consistency of the Zachary Karate Club network
#' sigma_c <- structural_consistency(g = karate_club, p_H = 0.1, epochs = 100, 
#'               k = round(sqrt(igraph::vcount(karate_club))))
#' avg_sigma_c <- mean(sigma_c)
#' 
#' @export
#' @importFrom igraph ecount vcount as_adjacency_matrix delete_edges E E<-
#' @import Matrix
#' @importFrom RSpectra eigs_sym
#' @importFrom dplyr tibble %>% arrange desc
#'
structural_consistency <- function(g, p_H = 0.1, epochs = 100, k = NA){
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  N <- vcount(g)
  L <- ecount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  
  # Allocate space for the result
  sigma_c <- numeric(epochs)
  
  for(i in 1:epochs){
    # Determine the link perturbation set dlt_E at random
    to_remove <- sample(L, links_to_remove)
    g_perturbed <- delete_edges(g, to_remove)
    
    # Compute the adjacency matrix of the graph after removal of dlt_E
    A_R <- as_adjacency_matrix(g_perturbed, type = "both")
    
    # Diagonalise A_R
    if(is.na(k) || k == N){
      XLX <- eigen(A_R, symmetric = T)
    }else{
      XLX <- eigs_sym(A_R, k = k, which = "LA")
    }
    
    # Compute the adjacency matrix of the perturbation set
    dlt_A <- A - A_R
    
    # Compute correction terms for eigenvalues based on dlt_A
    dlt_l <- rowSums((t(XLX$vectors) %*% dlt_A) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
    
    # Perturb A_R but keep eigenvectors unchanged
    A_perturbed <- XLX$vectors %*% diag(XLX$values + dlt_l) %*% t(XLX$vectors)
    
    non_edges <- get_non_edges(g_perturbed)
    
    prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2], 
                         spm = A_perturbed[non_edges]) %>% arrange(desc(spm))
    
    # Check if the top edges are indeed part of the original graph
    prediction <- prediction[1:links_to_remove, ]
    sigma_c[i] <- sum(g[from = prediction$nodeA, to = prediction$nodeB]) /
      links_to_remove
  }
  
  return(sigma_c)
}

