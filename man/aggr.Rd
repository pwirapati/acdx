\name{aggr}
\alias{aggr}
\title{Aggregate Single-Cell Expression Profiles}
\description{Produce aggregate matrix from raw counts}
\usage{
  aggr( x, cellgroups, v0 = 1/12, sep="\t", verbose=1 )
}
\arguments{
  \item{x}{Raw counts in sparse dgCMatrix format, with genes as rows
  and cells as columns. Alternatively, if it is a single string, it
  is interpreted as local file name containing text table (which
  can be compressed).}
  \item{cellgroups}{Labels for grouping the cells into the same aggregate,
    a vector with the same length as the number of columns of \code{x}}
  \item{v0}{Baseline variance, with the default 1/12 equals to that of
    the uniform distribution in an interval of width one}
  \item{sep}{Separator character; only if \code{x} is a file of text table}
  \item{verbose}{Output progress messages when performing the aggregation}
}
\details{
  The summary for each gene in an aggregate is the mean of raw counts,
  accompanied by an estimate of the variance of the mean (i.e., squared
  standard error). The aggregate data is a 3d array with the pair
  of mean and variance in the first dimension, the aggregate (as
  identified by \code{cellgroups} in the second, and the genes
  as the third.
}
\value{
  A list of:
  \item{annot}{Aggregate annotation which contains the aggregate
    identifier (\code{id}) and the number of cells per aggregate (\code{N}).}
  \item{y}{The aggregate array}
}
\seealso{
  \code{\link{acdx}}
}

\examples{
  \dontrun{
  ## downloads large files from GEO
  baseurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
  gsesuppl <- "GSE144nnn/GSE144735/suppl/"
  download.file(paste0(baseurl,gsesuppl,
    "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"),
    "KUL3.annotation.gz")
  download.file(paste0(baseurl,gsesuppl,
    "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"),
    "KUL3.rawUMI.gz")

  a <- read.delim("KUL3.annotation.gz")
  KUL3 <- aggr("KUL3.rawUMI.gz", paste(a$Cell_type,a$Sample,sep="|"))

  }

}
