# AC&#x26A1;DX: Aggregated Cells Differential eXpression

This is an R package for performing aggregation of single-cell expression profiles and rank genes based on specified patterns of changes in expression across combinations of cell types and conditions.

1. Take as input raw counts and assumed cell labels
2. Aggregate single cell profiles of the same cell label in the same sample, summarizing each gene by its mean expression of the raw counts and its standard error. The result is a 3D array of (mean,SE<sup>2</sup>) x aggregates x genes. The number of cells are also kept (common to all genes in an aggregate). 
3. Fit "meta-analytical" gamma GLM model which incorporate variance of the response values, included as offset in a variance component model that include gamma (multiplicative dispersion) terms. This is a joint GLM that simultaneously fit the mean (gene effects, sample effects) and dispersion (group-specific coefficients of variation, single-cell-level variance). Block bootstrapping by samples or patients is integral part of the implementation. The results are fitted model parameters x genes x bootstrap replicates.
4. Query the fitted models using various _interest functionals_. This is a generalization of _linear contrasts_ which can use non-linear functionals of the parameter space. A ranked gene list is produced, together with the significance statistics (FWER and FDP) that take into accounts multiple testing in the presence of dependence in genes and aggregates.
5. Visualize and interpret the ranked genes for a given query:
    - "top" table (similar to that of `limma`)
    - per-gene forest plots showing the aggregate summaries of top-ranking genes


[Tutorial](https://pwirapati.github.io/acdx/inst/doc/tutorial.html)
