
globalVariables(c("frac_rem", "method", "perf", "avg_perf", "links_rem", 
                  "links_recov"))

#' Performance evaluation of link predictors
#' 
#' Given a network of interest and any number of link prediction function names,
#' evaluates the performance of these link predictors by pruning edges from the 
#' network and checking if the methods give high likelihood scores to the
#' removed edges. Performance is assessed with four different metrics: Recall@k,
#' AUPR, AUROC and Average Precision.
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
#' \item{recall_at_k}{Fraction of top candidates links that are in the set of 
#' removed edges.}
#' \item{aupr}{The area under the Precision-Recall curve.}
#' \item{auroc}{The area under the Receiving Operating Characteristic curve.}
#' \item{avg_prec}{The average Precision at the point where Recall reaches its 
#' maximum value of 1.}
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
#' @importFrom purrr map2 map_dbl
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
                recall_at_k = numeric(col_len),
                aupr = numeric(col_len),
                auroc = numeric(col_len),
                avg_prec = numeric(col_len))
  
  evaluation <- map2(as.list(res$method), as.list(res$links_rem),
                         prune_predict_assess, g = g, 
                         removable_links = removable_links)
  
  res$recall_at_k <- map_dbl(evaluation, function(x) x["recall_at_k"])
  res$aupr <- map_dbl(evaluation, function(x) x["aupr"])
  res$auroc <- map_dbl(evaluation, function(x) x["auroc"])
  res$avg_prec <- map_dbl(evaluation, function(x) x["avg_prec"])
  
  return(res)
}

#' Determine the performance of a link predictor
#' 
#' Given a network of interest, the name of a link predictor, a number of links 
#' to prune from the network and a set of removable links, samples edges at 
#' random from the removable set and prunes the network. Then applies the link 
#' prediction method to the pruned topology and evaluates the prediction with 
#' four different metrics: Recall@k, AUPR, AUROC and Average Precision.
#' 
#' @param g igraph; A network of interest.
#' @param method character; The name of a link prediction function.
#' @param links_rem integer; Number of links to prune from the network.
#' @param removable_links igraph edge sequence; Set of edges that can be pruned.
#' 
#' @return A named numeric vector with the performance of the link predictor.
#' The elements of the vector are:
#' \item{recall_at_k}{Fraction of top candidates links that are in the set of 
#' removed edges.}
#' \item{aupr}{The area under the Precision-Recall curve.}
#' \item{auroc}{The area under the Receiving Operating Characteristic curve.}
#' \item{avg_prec}{The average Precision at the point where Recall reaches its 
#' maximum value of 1.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @importFrom igraph delete_edges
#' @import precrec
#' 
prune_predict_assess <- function(g, method, links_rem, removable_links){
  # Prepare vector for results
  perf <- numeric(4)
  names(perf) <- c("recall_at_k", "aupr", "auroc", "avg_prec")
  
  # Sample edges at random from the removable set
  edg <- sample(removable_links, links_rem)
  
  # Prune the network of interest
  g_perturbed <- delete_edges(g, edg)
  
  # Predict links and determine scores and classes
  pred <- do.call(method, list(g = g_perturbed))
  scr <- nrow(pred):1
  classes <- g[from = pred$nodeA, to = pred$nodeB]
  
  # Compute Recall@k
  perf["recall_at_k"] <- sum(classes[1:links_rem])/links_rem
  
  # Compute AUPR and AUROC
  meval <- evalmod(scores = scr, labels = classes)
  perf["aupr"] <- attr(meval$prcs[[1]], "auc")
  perf["auroc"] <- attr(meval$rocs[[1]], "auc")
  
  # Compute Avg. Precision
  perf["avg_prec"] <- mean(meval$prcs[[1]]$y[meval$prcs[[1]]$orig_points])
  
  return(perf)
}

#' Plot the performance evaluation of link predictors
#' 
#' Given the assessment of one or more link predictors with function 
#' \code{prune_recover}, plots the performance curves of the analysed methods
#' using the chose metric. Each point corresponds to the average performance for
#' a specific fraction of removed links across all considered epochs.
#' 
#' @param res tibble; The result of applying \code{prune_recover} to a network.
#' @param metric character; The metric that we want to plot. Should be one of 
#' 'recall_at_k', 'aupr', 'auroc' or 'avg_prec'.
#' @param colours character/numeric; A vector of colours to depict the 
#' performance curve of each link prediction method. There should be as many
#' colours as assessed link predictors.
#' @param err character; Include or not error bars in the plot. It should be one
#' of 'none', 'sd' or 'se'.
#' 
#' @return A ggplot with the performance curve(s) of the assesed prediction 
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
#' # Plot the performance curves of the considered link predictors
#' perf <- plot_lp_performance(assessment, colours = colours, err = "sd")
#' 
#' @export
#' @importFrom dplyr %>% mutate group_by summarise n
#' @import ggplot2
#' @importFrom stats sd
#'
plot_lp_performance <- function(res, metric = "recall_at_k", colours = NA, 
                                err = "none"){
  y_axis <- switch(metric, 
                   recall_at_k = "Recall@k", 
                   aupr = "AUPR", 
                   auroc = "AUROC", 
                   avg_prec = "Avg. Precision",
                   stop("The selected metric is not valid!")
                   )
  if(any(is.na(colours))){
    colours = seq_along(unique(res$method))
  }else if(length(colours) != length(unique(res$method))){
    stop("The number of link prediction methods and colours does not match!!!")
  }
  if(!(err %in% c("none", "sd", "se"))){
    stop("The indicated error type is not valid!")
  }
  
  res <- res %>% mutate(perf = res[[metric]])
  
  if(err == "none"){
    res <- res %>% group_by(method, frac_rem) %>% 
      summarise(avg_perf = mean(perf))
    err_bar_width <- 0
  }else if(err == "sd"){
    res <- res %>% group_by(method, frac_rem) %>% 
      summarise(avg_perf = mean(perf), 
                up_err = avg_perf + sd(perf), 
                dw_err = avg_perf - sd(perf))
    err_bar_width <- 0.01
  }else if(err == "se"){
    res <- res %>% group_by(method, frac_rem) %>% 
      summarise(avg_perf = mean(perf), 
                up_err = avg_perf + sd(perf)/sqrt(n()), 
                dw_err = avg_perf - sd(perf)/sqrt(n()))
    err_bar_width <- 0.01
  }
  
  p <- res %>% ggplot(aes_(~frac_rem, ~avg_perf, colour = ~method)) + 
    geom_line() + 
    geom_errorbar(aes_(ymin = ~dw_err, ymax = ~up_err), width = err_bar_width) +
    scale_colour_manual(values = colours) +
    labs(x = "Fraction of removed links", 
         y = y_axis, colour = "Predictor") + 
    theme_bw() + theme(legend.background = element_blank(), 
                       legend.position = "top")
  
  return(p)
}