
globalVariables(c("frac_rem", "method", "prec", "perf", "links_rem", 
                  "links_recov"))

#' Performance evaluation of link predictors
#' 
#' Given a network of interest and any number of link prediction function names,
#' assesses the performance of these link predictors by pruning edges from the 
#' network and checking if the methods give high likelihood scores to the
#' removed edges.
#' 
#' @param g igraph; The network of interest.
#' @param ... character; One or more function names able to produce a tibble
#' with 3 columns, the first two being node IDs (candidate link) and the third
#' their likelihood of interaction. This tibble should be sorted from most to 
#' less likely candidate links.
#' @param probes numeric; A vector indicating the fraction of links to be pruned
#' from \code{g} at random.
#' @param epochs integer; Number of times each fraction of links is randomly
#' removed.
#' @param preserve_conn logical; Should network connectivity be preserved after 
#' edge pruning? This is important for some prediction techniques that can only 
#' be applied to connected topologies. To preserve connectivity, the network's 
#' Minimum Spanning Tree (MST) is computed and links are pruned around it.
#' @param use_weights logical; If connectivity is to be preserved, this 
#' parameter indicates whether the MST should be computed taking edge weights 
#' into account.
#' 
#' @return Tibble with the following columns:
#' \item{method}{Link prediction method.}
#' \item{epoch}{Epoch of link removal.}
#' \item{frac_rem}{Fraction of links pruned at the specified epoch.}
#' \item{links_rem}{Number of links pruned at the specified epoch.}
#' \item{links_recov}{Number of removed links recoverd by the link predictor.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Lu, L. et al. (2015) Toward link predictability of complex 
#' networks. \emph{PNAS} 112(8):2325-2330
#'
#' @examples
#' # Assess the performance of three link predictors applied to the Zachary 
#' # Karate Club network
#' assessment <- prune_recover(g = karate_club, "lp_cn", "lp_aa", "lp_pa")
#' 
#' @export
#' @importFrom igraph mst E head_of tail_of
#' @importFrom dplyr tibble
#' @importFrom purrr map2_dbl
#'
prune_recover <- function(g, ..., probes = seq(0.1, 0.5, 0.1), epochs = 10,
                          preserve_conn = FALSE, use_weights = FALSE){
  predictors <- list(...)
  if(("lp_mce" %in% predictors || "lp_leig" %in% predictors || 
     "lp_isomap" %in% predictors) && !preserve_conn){
    stop(paste0("One or more of your link predictors requires a connected ",
                "network, set preserve_conn to TRUE!"))
  }
  if(preserve_conn){
    # Compute a minimum spanning tree and prune around it
    if(use_weights){
      g_mst <- mst(g)  
    }else{
      g_mst <- mst(g, algorithm = "unweighted") 
    }
    removable_links <- E(g)[-g[from = head_of(g_mst, E(g_mst)), 
                               to = tail_of(g_mst, E(g_mst)), edges = TRUE]]
  }else{
    removable_links <- E(g)
  }
  
  # Prepare tibble for the result
  col_len <- length(predictors) * length(probes) * epochs
  res <- tibble(method = rep(unlist(predictors), length(probes) * epochs), 
                epoch = rep(1:epochs, each = length(predictors)*length(probes)), 
                frac_rem = rep(rep(probes, each = length(predictors)), epochs), 
                links_rem = round(frac_rem * length(removable_links)), 
                links_recov = numeric(col_len))
  
  res$links_recov <- map2_dbl(as.list(res$method), as.list(res$links_rem),
                               prune_predict_assess, g = g, 
                               removable_links = removable_links)
  
  return(res)
}

#' Determine the performance of a link predictor
#' 
#' Given a network of interest, the name of a link predictor, a number of links 
#' to prune from the network and a set of removable links, samples edges at 
#' random from the removable set and prunes the network. Then applies the link 
#' prediction method to the pruned topology and checks if the 
#' top-\code{links_rem} ranked links are in the set of removed edges.
#' 
#' @param g igraph; A network of interest.
#' @param method character; The name of a link prediction function.
#' @param links_rem integer; Number of links to prune from the network.
#' @param removable_links igraph edge sequence; Set of edges that can be pruned.
#' 
#' @return The number of top-\code{links_rem} ranked links that are in the set
#' of removed edges.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @importFrom igraph delete_edges
#' 
prune_predict_assess <- function(g, method, links_rem, removable_links){
  # Sample edges at random from the removable set
  edg <- sample(removable_links, links_rem)
  
  # Prune the network of interest
  g_perturbed <- delete_edges(g, edg)
  
  # Predict links
  pred <- do.call(method, list(g = g_perturbed)) 
  pred <- pred[1:links_rem, ]
  
  # Determine the number of top predictions that are in the set of removed edges
  recovered <- sum(g[from = pred$nodeA, to = pred$nodeB])
  
  return(recovered)
}

#' Plot the performance evaluation of link predictors
#' 
#' Given the assessment of one or more link predictors with function 
#' \code{prune_recover}, plots the precision curves of the analysed methods.
#' Each point corresponds to the average precision for a specific fraction of
#' removed links across all considered epochs.
#' 
#' @param res tibble; The result of applying \code{prune_recover} to a network.
#' @param colours character/numeric; A vector of colours to depict the 
#' performance curve of each link prediction method. There should be as many
#' colours as assessed link predictors.
#' 
#' @return A ggplot with the precision curve(s) of the assesed prediction 
#' method(s).
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Assess the performance of three link predictors applied to the Zachary 
#' # Karate Club network
#' assessment <- prune_recover(g = karate_club, "lp_cn", "lp_aa", "lp_pa")
#' 
#' # Define a set of colours to plot
#' colours <- c("#8da0cb", "#fc8d62", "#66c2a5")
#' 
#' # Plot the precision curves of the considered link predictors
#' perf <- plot_lp_precision(assessment, colours)
#' 
#' @export
#' @importFrom dplyr %>% mutate group_by summarise
#' @import ggplot2
#'
plot_lp_precision <- function(res, colours = NA){
  if(any(is.na(colours))){
    colours = seq_along(unique(res$method))
  }else if(length(colours) != length(unique(res$method))){
    stop("The number of link prediction methods and colours does not match!!!")
  }
  
  p <- res %>% mutate(prec = links_recov/links_rem) %>% 
    group_by(method, frac_rem) %>% 
    summarise(perf = mean(prec)) %>% 
    ggplot(aes_(~frac_rem, ~perf, colour = ~method)) + 
    geom_line() + 
    scale_colour_manual(values = colours) +
    labs(x = "Fraction of removed links", 
         y = "Precision", colour = "Predictor") + 
    theme_bw() + theme(legend.background = element_blank(), 
                       legend.position = "top")
  
  return(p)
}

#' Plot the performance evaluation of link predictors
#' 
#' Given the assessment of one or more link predictors with function 
#' \code{prune_recover}, plots the precision curves of the analysed methods.
#' Each point corresponds to the average precision for a specific fraction of
#' removed links across all considered epochs. This function also plots error 
#' bars corresponding to standard errors.
#' 
#' @param res tibble; The result of applying \code{prune_recover} to a network.
#' @param colours character/numeric; A vector of colours to depict the 
#' performance curve of each link prediction method. There should be as many
#' colours as assessed link predictors.
#' 
#' @return A ggplot with the precision curve(s) of the assesed prediction 
#' method(s).
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Assess the performance of three link predictors applied to the Zachary 
#' # Karate Club network
#' assessment <- prune_recover(g = karate_club, "lp_cn", "lp_aa", "lp_pa")
#' 
#' # Define a set of colours to plot
#' colours <- c("#8da0cb", "#fc8d62", "#66c2a5")
#' 
#' # Plot the precision curves of the considered link predictors
#' perf <- plot_lp_precision_se(assessment, colours)
#' 
#' @export
#' @importFrom dplyr %>% mutate group_by summarise n
#' @import ggplot2
#' @importFrom stats sd
#'
plot_lp_precision_se <- function(res, colours = NA){
  if(any(is.na(colours))){
    colours = seq_along(unique(res$method))
  }else if(length(colours) != length(unique(res$method))){
    stop("The number of link prediction methods and colours does not match!!!")
  }
  
  p <- res %>% mutate(prec = links_recov/links_rem) %>% 
    group_by(method, frac_rem) %>% 
    summarise(perf = mean(prec), 
              up_err = perf + sd(prec)/sqrt(n()), 
              dw_err = perf - sd(prec)/sqrt(n())) %>% 
    ggplot(aes_(~frac_rem, ~perf, colour = ~method)) + 
    geom_line() + 
    geom_errorbar(aes_(ymin = ~dw_err, ymax = ~up_err), width = 0.01) +
    scale_colour_manual(values = colours) +
    labs(x = "Fraction of removed links", 
         y = "Precision", colour = "Predictor") + 
    theme_bw() + theme(legend.background = element_blank(), 
                       legend.position = "top")
  
  return(p)
}

#' Plot the performance evaluation of link predictors
#' 
#' Given the assessment of one or more link predictors with function 
#' \code{prune_recover}, plots the precision curves of the analysed methods.
#' Each point corresponds to the average precision for a specific fraction of
#' removed links across all considered epochs. This function also plots error 
#' bars corresponding to standard deviations.
#' 
#' @param res tibble; The result of applying \code{prune_recover} to a network.
#' @param colours character/numeric; A vector of colours to depict the 
#' performance curve of each link prediction method. There should be as many
#' colours as assessed link predictors.
#' 
#' @return A ggplot with the precision curve(s) of the assesed prediction 
#' method(s).
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Assess the performance of three link predictors applied to the Zachary 
#' # Karate Club network
#' assessment <- prune_recover(g = karate_club, "lp_cn", "lp_aa", "lp_pa")
#' 
#' # Define a set of colours to plot
#' colours <- c("#8da0cb", "#fc8d62", "#66c2a5")
#' 
#' # Plot the precision curves of the considered link predictors
#' perf <- plot_lp_precision_sd(assessment, colours)
#' 
#' @export
#' @importFrom dplyr %>% mutate group_by summarise
#' @import ggplot2
#' @importFrom stats sd
#'
plot_lp_precision_sd <- function(res, colours = NA){
  if(any(is.na(colours))){
    colours = seq_along(unique(res$method))
  }else if(length(colours) != length(unique(res$method))){
    stop("The number of link prediction methods and colours does not match!!!")
  }
  
  p <- res %>% mutate(prec = links_recov/links_rem) %>% 
    group_by(method, frac_rem) %>% 
    summarise(perf = mean(prec), 
              up_err = perf + sd(prec), 
              dw_err = perf - sd(prec)) %>% 
    ggplot(aes_(~frac_rem, ~perf, colour = ~method)) + 
    geom_line() + 
    geom_errorbar(aes_(ymin = ~dw_err, ymax = ~up_err), width = 0.01) +
    scale_colour_manual(values = colours) +
    labs(x = "Fraction of removed links", 
         y = "Precision", colour = "Predictor") + 
    theme_bw() + theme(legend.background = element_blank(), 
                       legend.position = "top")
  
  return(p)
}