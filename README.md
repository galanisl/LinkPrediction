LinkPrediction
================
Gregorio Alanis-Lobato

Implementation of different classes of link predictors and methods to assess and visualise their performance.

Introduction
============

Complex systems can be represented by graphs *G*(*V*, *E*), where *E* is the set of interactions (edges) between a set *V* of system components (nodes). Link predictors assign likelihood scores of interaction to all node pairs that are disconnected in the observable network topology. These scores can lead to the prediction of future friendships in a social network or future direct flights between cities in an airport transportation network.

To assess the performance of link predictors, one normally removes an increasing number of edges *E*<sup>*R*</sup> ⊂ *E* from the network of interest, applies a link prediction method to the pruned network and checks whether the |*E*<sup>*R*</sup>| candidate links with highest likelihood scores are in *E*<sup>*R*</sup>.

`LinkPrediction` implements the following classes of link prediction techniques:

-   Neighbourhood-based predictors, which assign high likelihood scores to node pairs that share many neighbours (`lp_cn`, `lp_pa`, `lp_jc`, `lp_dice`, `lp_aa` and `lp_ra`).
-   CAR-based predictors, which assign high likelihood scores to node pairs that share many neighbours that also interact with themselves (`lp_car`, `lp_cpa`, `lp_caa` and `lp_cra`).
-   Embedding-based methods rank non-adjacent nodes through distances on a network projection in two-dimensional Euclidean space (`lp_isomap`, `lp_leig`, and `lp_mce`).
-   The Hierarchical Random Graph method, which searches the space of all possible dendrograms of a network for the ones that best fit its hierarchical structure. Non-adjacent node pairs that have high average probability of being connected within these dendrograms are considered good candidates for interaction (`lp_hrg`).
-   The Structural Perturbation Method, which is based on the hypothesis that links are predictable if removing them has only small effects on network structure (`lp_spm`).

In addition, `LinkPrediction` includes function `lp_matrix` and `get_non_edges`, which are key to extend the package with other link prediction approaches (more on this below).

Finally, the performance of these link predictors on a network of interest can be assessed and visualised with functions `prune_recover` and `plot_lp_precision`.

For more information on the link prediction problem, see:

-   Martínez, V., Berzal, F. & Cubero, J.-C. A survey of link prediction in complex networks. *ACM Computing Surveys* **49**, 1-33 (2016).
-   Lü, L., Pan, L., Zhou, T., Zhang, Y.-C. & Stanley, H. E. Toward link predictability of complex networks. *PNAS* **112**, 2325-2330 (2015).

Installation
============

1.  Install the `devtools` package from CRAN if you haven't done so:

``` r
install.packages("devtools")
```

1.  Load the `devtools` package:

``` r
library("devtools")
```

1.  Install `LinkPrediction` using the `install_github` function:

``` r
install_github("galanisl/LinkPrediction")
```

Usage
=====

To start using `LinkPrediction`, load the package:

``` r
library("LinkPrediction")
```

`LinkPrediction` includes two complex network datasets (type `?karate_club` and `?jazz_collab` in R for more information). Let's use the Jazz Collaboration network to illustrate the use of the package.

First, let's apply the Common Neighbours link predictor to the network:

``` r
(cn <- lp_cn(jazz_collab))
```

    ## # A tibble: 16,761 x 3
    ##    nodeA nodeB   scr
    ##    <int> <int> <dbl>
    ##  1     7   111    41
    ##  2    27    79    37
    ##  3    27    66    33
    ##  4     7   113    29
    ##  5    66   118    29
    ##  6    74    89    28
    ##  7    88    92    27
    ##  8    23   124    27
    ##  9    53    66    26
    ## 10    27    99    26
    ## # ... with 16,751 more rows

Note how the lp\_\* functions return a `tibble` of candidate links with their likelihoods of interaction. Let's now compute the structural consistency *σ*<sub>*c*</sub> ∈ \[0, 1\] of this network. The higher it is, the higher its link predictability:

``` r
sigma_c <- structural_consistency(jazz_collab)
mean(sigma_c)
```

    ## [1] 0.705438

``` r
sd(sigma_c)
```

    ## [1] 0.01962508

Since the *σ*<sub>*c*</sub> is high, link predictors are likely to give us good candidates of interaction.

Let's now implement a link predictor of our own with the help of function `get_non_edges`. Our method will return a tibble with a random ordering of the disconnected node pairs in the network:

``` r
lp_rnd <- function(g){
  non_edges <- get_non_edges(g)
  prediction <- tibble(nodeA = non_edges[, 1], nodeB = non_edges[, 2],
                       scr = sample(nrow(non_edges))) %>% 
    arrange(desc(scr))
  return(prediction)
}
```

Finally, let's compare our method with other link predictors and visualise the results:

``` r
assessment <- prune_recover(jazz_collab, 
                            "lp_rnd", 
                            "lp_cn", 
                            "lp_car", 
                            "lp_spm")

plot_lp_precision_sd(assessment)
```

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png)
