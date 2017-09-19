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

    ## [1] 0.7076642

``` r
sd(sigma_c)
```

    ## [1] 0.01974268

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

Session information
===================

``` r
sessionInfo()
```

    ## R version 3.4.1 (2017-06-30)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bindrcpp_0.2       LinkPrediction_0.9 dplyr_0.7.3       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] igraph_1.1.2     Rcpp_0.12.12     knitr_1.16       bindr_0.1       
    ##  [5] magrittr_1.5     munsell_0.4.3    colorspace_1.3-2 lattice_0.20-35 
    ##  [9] R6_2.2.2         rlang_0.1.2      plyr_1.8.4       stringr_1.2.0   
    ## [13] tools_3.4.1      grid_3.4.1       gtable_0.2.0     RSpectra_0.12-0 
    ## [17] htmltools_0.3.6  lazyeval_0.2.0   yaml_2.1.14      assertthat_0.2.0
    ## [21] rprojroot_1.2    digest_0.6.12    tibble_1.3.4     Matrix_1.2-10   
    ## [25] purrr_0.2.3      ggplot2_2.2.1    glue_1.1.1       evaluate_0.10.1 
    ## [29] rmarkdown_1.6    labeling_0.3     stringi_1.1.5    compiler_3.4.1  
    ## [33] scales_0.4.1     backports_1.1.0  pkgconfig_2.0.1
