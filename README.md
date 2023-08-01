# xcore

xcore is an R package for transcription factor activity modeling
based on known molecular signatures and user's gene expression data.
Accompanying [xcoredata](https://github.com/mcjmigdal/xcoredata/) package
provides a collection of molecular signatures, constructed from publicly
available ChiP-seq experiments.

We refer interested users to our publication in [BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05084-0).

## Installation
xcore and xcoredata can be installed from Bioconductor:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("xcore")
BiocManager::install("xcoredata")
```

## Usage

A vignette showing xcore basic usage is available [here](https://bkaczkowski.github.io/xcore/articles/xcore_vignette.html).

## Parallel computing

xcore can take advantage of parallelization to speed up calculations, especially for model
training and estimates testing. To use parallel computing in `R` one have to first
register parallel backend. While there are many parallel backends to choose 
from, internally xcore uses [`foreach`](https://cran.r-project.org/web/packages/foreach)
to implement parallel computing. Having this in mind we should use a backend
supported by `foreach`. 

In the vignette we are using [`doParallel`](https://cran.r-project.org/package=doParallel)
backend, together with [`BiocParallel`](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
package providing unified interface across different OS. Those packages can be
installed with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocParallel")
install.packages("doParallel")
```
