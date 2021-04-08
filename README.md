# rlingam: R implementation of LiNGAM algorithms

## Implemented Algorithms
- ICALiNGAM
- DirectLiNGAM

## Install
```r
remotes::install_github("gkikuchi/rlingam")
```

## Usage
```r
library(rlingam)

X <- gen_dummy_data(random_state = 10)

# icalingam
mdl <- ICALiNGAM$new()
mdl$fit(X)

# directlingam
mdl <- DirectLiNGAM$new()
mdl$fit(X)

print(mdl$causal_order)
print(mdl$adjacency_matrix)
plot_adjacency_mat(mdl$adjacency_matrix, node_labels = names(X))
```
