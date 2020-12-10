# AC&#x26A1;DX: Aggregated Cells Differential eXpression

## NOTE

There is major change of interface in v1.1.0. If you'd like to keep old analysis, install by specifying the version:
```
devtools::install_github("pwirapati/acdx@v1.1.0")
```


## Blurbs

This is an R package for performing aggregation of single-cell expression profiles and rank genes based on specified patterns of changes in expression across combinations of cell types and conditions.

Feature highlights:

* Aggregation of cell expression profiles based on average and standard error of the raw counts (instead of pseudo-bulk sums)
* Random-effect meta-regression with gamma-GLM-like error model, with fast implementation to be easily applied to all genes, all cell types and many bootstrap replicates
* Ranking based on flexible _interest functionals_ that can span across cell types and conditions, allowing selecting genes based on interesing biological relationship across cell types
* Multiple-testing significance analysis based on bootstrap resampling
* Visualization of per-gene expression patterns using forest plots

[Tutorial for v1.1.0](https://pwirapati.github.io/acdx/inst/doc/tutorial.html)
