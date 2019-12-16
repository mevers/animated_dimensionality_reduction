# animated_dimensionality_reduction

Animation showing the results of three different dimensionality reduction (DR) techniques in R.

## Overview

We compare results from three different DR techniques

- principle component analysis (PCA),
- t-distributed stochastic neighbour embedding (t-SNE), and
- uniform manifold approximation and projection (UMAP)

We don't give any details involving the different DR techniques, as this would go beyond the scope of this short example and wouldn't do justice to the plethora of excellent material available online. Some references with links are given at the end of this README.md.

As source data we use a subset of 10,000 images from the [MNIST database of handwritten digits](http://yann.lecun.com/exdb/mnist/). The original MNIST training database contains 60,000 28x28 pixel images from approximately 250 writers. The input data consists of a 10,000x784 matrix.

We use `tweenr` to create a smooth animation linking results from the three different DR techniques. [`tweenr`](https://github.com/thomasp85/tweenr) creates smooth animations by calculating intermediate values between different states (here, the different DR results). We then use `animation::saveGIF` to store the animation in an animated GIF.


## Prerequisities

In order to use t-SNE and UMAP in R the packages `Rtsne` and `umap` must be installed:

```r
install.packages("Rtsne")
install.packages("umap")
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
df <- read_csv("train.zip")
df <- df %>% head(10000)
```

The resulting `df` is a 10,000x785-dimensional `data.frame`, where the first column `df[, 1]` contains the label (the number) and the remaining 784 columns encode the 28x28 pixels.

## Perform dimensionality reduction

We perform dimensionality reduction using the functions

- `prcomp` for PCA,
- `Rtsne::Rtsne` for t-SNE, and
- `umap::umap` for UMAP.

All functions take the 10,000x784 matrix `df[, -1]` as input. We store results in a `res_lst` and timing information from `system.time` in `t_lst`.

```r
library(Rtsne)
library(umap)
t_lst <- list()
res_lst <- list()
t_lst[["t-SNE"]] <- system.time(
    res_lst[["t-SNE"]] <- Rtsne::Rtsne(
        df[, -1], dims = 2, perplexity = 30, max_iter = 500))
t_lst[["UMAP"]] <- system.time(
    res_lst[["UMAP"]] <- umap::umap(df[, -1]))
t_lst[["PCA"]] <- system.time(
    res_lst[["PCA"]] <- prcomp(df[, -1]))
```

## Timings

The timing results of the three different DR techniques are shown in the following figure.

```r
t_df <- t_lst %>%
    map(stack) %>%
    bind_rows(.id = "algorithm")

t_df %>%
    pivot_wider(values_from = values, names_from = ind)
## A tibble: 3 x 6
#  algorithm user.self sys.self elapsed user.child sys.child
#  <chr>         <dbl>    <dbl>   <dbl>      <dbl>     <dbl>
#1 t-SNE          68.4    1.57     85.2          0         0
#2 UMAP          104.    10.0     132.           0         0
#3 PCA            25.4    0.464    28.2          0         0    

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
    keep_state(10) %>%
    tween_state(
        filter(data, algorithm == "t-SNE"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(10) %>%
    tween_state(
        filter(data, algorithm == "UMAP"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(10) %>%
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
