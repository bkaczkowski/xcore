# 4/29/22
1.1.1 xcore is now released on Bioconductor

# 8/24/22
1.2.0 the following changes and additions has been made:
+ using the replicate average and replicate pooled standard errors has been
  changed to using weighted average with variance defined weights and the 
  standard error of weighted average. The use of latter has better theoretical
  grounding and its use is evidenced in literature (Arner E, et al. PMC3402332).
+ a bug causing each of the models produced by `modelGeneExpression` to carry
  a hard copy of molecular signatures matrix (X) and leading to large objects
  sizes has been resolved.
+ `translateCounts` function has been added to ease construction of gene level
  analysis.
+ `regressionData` function has been added allowing to construct analyses using
  customly pre-processed data.
+ new vignette has been added showing how to conduct a gene level analysis.