# animated_dimensionality_reduction

Animation showing the results of different dimensionality reduction (DR) techniques in R.


## Overview

We compare results from four different DR techniques

- principle component analysis (PCA),
- t-distributed stochastic neighbour embedding (t-SNE),
- uniform manifold approximation and projection (UMAP), and
- LargeVis

We don't give any details involving the different DR techniques, as this would go beyond the scope of this short example and wouldn't do justice to the plethora of excellent material available online. Some references with links are given at the end of this [README](README.md).

As source data we use a subset of 10,000 images from the [MNIST database of handwritten digits](http://yann.lecun.com/exdb/mnist/). The original MNIST training database contains 60,000 28x28 pixel images from approximately 250 writers. Here the input data consist of a 10,000x784 matrix.

We use `tweenr` to create a smooth animation linking results from the three different DR techniques. [`tweenr`](https://github.com/thomasp85/tweenr) creates smooth animations by calculating intermediate values between different states (here, the different DR results). We then use `animation::saveGIF` to store the animation in an animated GIF.


## Prerequisities

In order to use t-SNE and UMAP in R the packages `Rtsne` and `umap` must be installed (available from CRAN); `largeVis` is [available from GitHub](https://github.com/elbamos/largeVis) and can be installed via `devtools::install_github`.

```r
install.packages("Rtsne")
install.packages("umap")
devtools::install_github("elbamos/largeVis")
```

Furthermore we need the following packages for processing data and writing the final animated GIF.

```r
install.packages("tidyverse")
install.packages("hrbrthemes")
install.packages("animation")
```


## Loading the data

We load the MNIST traning data (which has already been [conveniently formatted as a CSV file](https://pjreddie.com/projects/mnist-in-csv/)), and select the first 10,000 rows corresponding to the first 10,0000 images from the MNIST training database.

```r
library(tidyverse)
df <- read_csv("train.zip") %>% head(10000)
```

The resulting `df` is a 10,000x785-dimensional `data.frame`, where the first column `df[, 1]` contains the label (the number) and the remaining 784 columns encode the 28x28 pixels.

## Perform dimensionality reduction

We perform dimensionality reduction using the functions

- `prcomp` for PCA,
- `Rtsne::Rtsne` for t-SNE,
- `umap::umap` for UMAP, and
- `largeVis::largeVis` for largeVis

All functions take the 10,000x784 matrix `df[, -1]` as input. We store results in a `res_lst` and timing information from `system.time` in `t_lst`.

```r
library(Rtsne)
library(umap)
library(largeVis)
t_lst <- list()
res_lst <- list()
t_lst[["t-SNE"]] <- system.time(
    res_lst[["t-SNE"]] <- Rtsne::Rtsne(
        df[, -1], dims = 2, perplexity = 30, max_iter = 500))
t_lst[["UMAP"]] <- system.time(
    res_lst[["UMAP"]] <- umap::umap(df[, -1]))
t_lst[["PCA"]] <- system.time(
    res_lst[["PCA"]] <- prcomp(df[, -1]))
t_lst[["largeVis"]] <- system.time(
    res_lst[["largeVis"]] <- largeVis(
        df[, -1], n_trees = 50, tree_th = 200, K = 50))
```


## Timings

The timing results of the different DR techniques are shown in the following figure.

```r
t_df <- t_lst %>%
    map(stack) %>%
    bind_rows(.id = "algorithm")

t_df %>%
    pivot_wider(values_from = values, names_from = ind)
## A tibble: 4 x 6
#  algorithm user.self sys.self elapsed user.child sys.child
#  <chr>         <dbl>    <dbl>   <dbl>      <dbl>     <dbl>
#1 t-SNE          36.0    0.791    37.5          0         0
#2 UMAP           58.0   12.1      70.2          0         0
#3 PCA            18.7    0.152    18.9          0         0
#4 largeVis      394.     0.809   394.           0         0    

library(hrbrthemes)
t_df %>%
    filter(values > 0) %>%
    rename(times = ind) %>%
    ggplot(aes(algorithm, values, fill = times)) +
    geom_col() +
    theme_ft_rc()
```

![Figure `timings.png` not found](./timings.png)


## Animation

We first define a helper function that extracts data of the first two reduced dimensions for every image and for every DR technique.  

```r
get_data <- function(res, labels) {
    if (class(res) == "umap") {
        df <- res$layout
    } else if (class(res) == "prcomp") {
        df <- res$x[, 1:2]
    } else if (class(res) == "largeVis") {
        df <- t(res$coords)
    } else if ("Y" %in% names(res)) df <- res$Y
    df %>%
        as.data.frame() %>%
        bind_cols(cbind.data.frame(Label = labels)) %>%
        setNames(c("Dimension 1", "Dimension 2", "Label")) %>%
        mutate(Label = as.factor(Label))
}
```

We then extract data from `res_lst`

```r
data <- res_lst %>%
    map(get_data, df[, 1]) %>%
    bind_rows(.id = "algorithm")
```

and create the "tweened" data

```r
library(tweenr)
data_tween <- filter(data, algorithm == "PCA") %>%
    keep_state(20) %>%
    tween_state(
        filter(data, algorithm == "t-SNE"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    tween_state(
        filter(data, algorithm == "UMAP"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    tween_state(
        filter(data, algorithm == "largeVis"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    tween_state(
        filter(data, algorithm == "PCA"),
        ease = "cubic-in-out", nframes = 100) %>%
    group_split(.frame)
```

`data_tween` is a `list` of `data.frame`s, one for every frame.

We now create the base `ggplot2`-based `grob`

```r
gg_base <- ggplot() +
    geom_text(
        aes(x = `Dimension 1`, y = `Dimension 2`, colour = Label, label = Label),
        size = 3, alpha = 0.5, show.legend = F) +
    theme_ft_rc()
```

We are now ready to generate the animated GIF by looping through the tweened data and adding the per-frame `data.frame`s.

```r
library(animation)
oopt <- ani.options(interval = 1 / 20)
i <- 1
saveGIF({
    for (d in data_tween) {
        cat(sprintf("Processing %i/%i\n", i, length(data_tween)))
        gg <- gg_base %+% d
        gg <- gg + labs(
            title = sprintf("Dimensionality reduction algorithm: %s",
            unique(d$algorithm)),
            subtitle = "Source data: Subset of MNIST",
            caption = "Author: Maurits Evers (maurits.evers@gmail.com)")
        plot(gg)
        i <- i + 1
    }},
    movie.name = "animation.gif",
    ani.width = 800, ani.height = 640)
```

![Figure `animation.gif` not found](./animation.gif)


## References

### t-SNE
- [L. van der Maaten and G. Hinton, *Visualising Data using t-SNE*, Journal of Machine Learning 9, 2579 (2008)  ](http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf)
- [Comprehensive Guide on t-SNE algorithm with implementation in R & Python](https://www.analyticsvidhya.com/blog/2017/01/t-sne-implementation-r-python/)
- [Quick and easy t-SNE analysis in R](https://www.r-bloggers.com/quick-and-easy-t-sne-analysis-in-r/)

### UMAP
- [L. McInnes, J. Healy and J. Melville, *UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction*, arXiv:1802.03426 (2018)](https://arxiv.org/abs/1802.03426)
- [UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction](https://umap-learn.readthedocs.io/en/latest/index.html)
- [Uniform Manifold Approximation and Projection in R](https://cran.r-project.org/web/packages/umap/vignettes/umap.html)

### largeVis
- [J. Tang et al., *Visualizing Large-scale and High-dimensional Data*, arXiv:1602.00370 (2016)](https://arxiv.org/abs/1602.00370)


## Comments

### A comment on the timings results

- `Rtsne` is a wrapper around the [fast C++ based Barnes-Hut implementation of the t-SNE algorithm](https://github.com/lvdmaaten/bhtsne/) from the original author (Laurens van der Maaten).
- By default, `umap` uses an implementation of the UMAP algorithm written in R. It is possible to use a faster Python-based implementation via the `umap-learn` Python package. Details on how to do that can be found on [Interfacing with python package ‘umap-learn’](https://cran.r-project.org/web/packages/umap/vignettes/umap_learn.html).


## `sessionInfo`

```r
sessionInfo()
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Sierra 10.12.6
#
#Matrix products: default
#BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
#
#locale:
#[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
# [1] animation_2.6    tweenr_1.0.1     hrbrthemes_0.7.3 umap_0.2.4.0
# [5] Rtsne_0.15       forcats_0.4.0    stringr_1.4.0    dplyr_0.8.3
# [9] purrr_0.3.3      readr_1.3.1      tidyr_1.0.0      tibble_2.1.3
#[13] ggplot2_3.2.1    tidyverse_1.2.1
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2        lubridate_1.7.4   lattice_0.20-38   assertthat_0.2.1
# [5] zeallot_0.1.0     digest_0.6.22     utf8_1.1.4        RSpectra_0.16-0
# [9] R6_2.4.0          cellranger_1.1.0  backports_1.1.5   evaluate_0.14
#[13] httr_1.4.1        pillar_1.4.2      gdtools_0.2.1     rlang_0.4.1
#[17] lazyeval_0.2.2    readxl_1.3.1      rstudioapi_0.10   extrafontdb_1.0
#[21] magick_2.2        Matrix_1.2-17     reticulate_1.13   rmarkdown_1.16
#[25] labeling_0.3      extrafont_0.17    munsell_0.5.0     broom_0.5.2
#[29] compiler_3.6.1    modelr_0.1.5      xfun_0.10         pkgconfig_2.0.3
#[33] askpass_1.1       systemfonts_0.1.1 htmltools_0.4.0   openssl_1.4.1
#[37] tidyselect_0.2.5  fansi_0.4.0       crayon_1.3.4      withr_2.1.2
#[41] grid_3.6.1        nlme_3.1-141      jsonlite_1.6      Rttf2pt1_1.3.7
#[45] gtable_0.3.0      lifecycle_0.1.0   magrittr_1.5      scales_1.0.0
#[49] cli_1.1.0         stringi_1.4.3     farver_2.0.1      xml2_1.2.2
#[53] generics_0.0.2    vctrs_0.2.0       tools_3.6.1       glue_1.3.1
#[57] hms_0.5.1         colorspace_1.4-1  rvest_0.3.4       knitr_1.25
#[61] haven_2.1.1
```
