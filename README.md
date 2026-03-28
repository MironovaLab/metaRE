# metaRE
R package for motif discovery via meta-analysis of microarray and RNA-seq data.

## System requirements

- R >= 3.3.2
- Compiler with C++14 support

## Installation

1. Install the CRAN and Bioconductor dependencies:
```r
install.packages(c("Rcpp", "BH", "futile.logger", "foreach"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("limma", "edgeR", "GEOquery", "Biobase"))
```

2. Download `metaRE`:
```
git clone https://github.com/cheburechko/metaRE
```

3. Install `metaRE` from the command line:
```
R CMD INSTALL [path/to/package]
```

Or from R using `devtools`:

```r
install.packages("devtools")
devtools::install("[path/to/package]")
```
