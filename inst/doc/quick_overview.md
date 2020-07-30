# Quick Overview

## Installation

At the moment, `acdx` is not in CRAN yet. 

> It will not be in Bioconductor because it does not depend on other Bioconductor packages.

Only the source package is available through GitHub, and it requires C compiler with OpenMP version 4.5 or above.
This should be available through the most recent version of Xcode (for Mac) or RTools (for Window).

> TODO: test and confirm by people with such setups

### Plain download or `git clone`

The zip archive of the source package can be downloaded from GitHub (either through the GUI or
by cloning the repo using git command line). R Package installation can then proceed
in the usual way (`R CMD ...` in shell or `install.packages()` inside R).

### Using `devtools`

If the package `devtools` is installed, it can be downloaded by

```
install_github("pwirapati/acdx")
```

## Data Preparation

Two cell-level data objects are required:

1. Raw read counts, can be ordinary R matrix or a sparse `dgCMatrix` from the `Matrix` package
2. Cell labels (a vector of `character` or `factor`), in the same order as the count matrix columns

As an example, we will use CRC data from GSE144735. Both data objects are available from GEO, in
the form of compressed tab-delimited (non-sparse) tables.

1. [Cell annotation](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144735&format=file&file=GSE144735%5Fprocessed%5FKUL3%5FCRC%5F10X%5Fannotation%2Etxt%2Egz) (190kb)
2. [Raw counts](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144735&format=file&file=GSE144735%5Fprocessed%5FKUL3%5FCRC%5F10X%5Fraw%5FUMI%5Fcount%5Fmatrix%2Etxt%2Egz)(54 Mb)

The annotation can be read in as data frame directly. The raw counts, if read in as data.frame, need to be converted into a matrix or sparse matrix.

## Aggregation

Load the package

```
library(acdx)
```

Assume that `annot` and `rumi` are cell-level annotation and raw UMI data.

Decide on the aggregation keys (cell groups).

This example aggregate into cell types, regardless of the samples
```
KUL3_ctype <- aggr( rumi, annot$Cell_type )
```

This example aggregate into a combination (Cartesian product) of cell types and sample id

```
KUL3_ctype_sample <- aggr( rumi, paste(annot$Cell_type,annot$Sample,sep="|") )
```

Let's examine the objects.

## Miscellanea

Aggregate variance:

<div>
\[
v_{ij} = \left[ \frac{1}{12} + s_{ij}^2 \right]/n_i
\]
</div>
