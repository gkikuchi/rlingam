# rlingam: R implementation of LiNGAM algorithms
Under Development

## Implemented Algorithms
- Direct LiNGAM

## Not Yet Implemented
- Other algorithms
- Bootstrapping
- etc

## Install
```r
remotes::install_github("gkikuchi/rlingam")
```

## Usage
```r
library(rlingam)

X <- rlingam::gen_dummy_data(random_state = 10)
mdl <- rlingam::DirectLiNGAM$new()
mdl$fit(X)

print(mdl$causal_order)
print(mdl$adjacency_matrix)
plot_adjacency_mat(mdl$adjacency_matrix, node_labels = names(X))
```
