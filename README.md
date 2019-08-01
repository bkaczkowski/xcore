# xcore

<!-- badges: start -->
<!-- badges: end -->

The goal of xcore is to ...

## To implement:

* reads to counts DONE
    - CTSS to counts module using data.table
    - BAM    to 5' counts module using Bioconductor packages
* annotation and collapsing of genomic features DONE
    - importing annotation (prepared by ChIPseeker package and Gencode/USCS knownGene annotations)
    - ENTREZID (gene) level collapsing 
    - F5 promoter and Enhancer level (genomic ranges level)
* Differential expression DONE/TESTING
    - DEseq DE
    - edgeR DE 
    - (optional) limma for MARA or other normally distributed scores/values
* GSEA worklow IN PROGRESS
    - ENTREZID (gene) level expression as input
    - DE using edgeR or DeSeq
    - GSEA using fgsea package
    - visualization of densely overlapping genes and genesets by pheatmap (similar to enrichplot::heatplot)
    - visualization of KEGG pathway using pathview::pathview()
* MARA workflow
    - export the count table and ranges (F5 promoter level) to MARA input format
    - write a wrapper to call MARA scripts
    - read in the MARA output and plot the most variant TFs
    - do DE analysis on MARA scores using limma and plot to differentially active TFs
* Chia-PET workflow

    
    
    

## Installation

You can install the released version of xcore from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("xcore")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(xcore)
## basic example code
```

