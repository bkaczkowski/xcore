#' edgeR wrapper
#'
#' @param counts data.matrix with the counts data, rownames should represent features and colnames sample names
#' @param sample_annot data.frame with sample annotation sample_annot$sample_names should be same as colnames(counts), it has to contain exp_fac column with experimental factor, control needs to be a the first level of the factor (see relevel() )
#' @param design design matrix
#' @param contrasts matrix with contrasts
#' @param robust should robust statitsital modeling be used?
#' @param plots should plotBCV() and plotQLDisp() be used to generate plots?
#' @param adj_pv_threshold adjusted p-value significance threshold for DE testing, 0.05 by default
#' @param log2fc_threshold log2 Fold Change threshold for DE testing, 1 by default#'
#' @return a list of differential expression results, each element of the list corresponds to one row of contrast matrix
#' @export
#' @importFrom  edgeR DGEList calcNormFactors estimateDisp glmQLFit plotBCV plotQLDisp
edgeR_DE_wrapper = function(counts, sample_annot, design, contrasts, robust = TRUE, plots = TRUE,
                            log2fc_threshold = 1,adj_pv_threshold = 0.05 ){

  #initial checks
  if ( sum ( sample_annot$sample_names != colnames(counts) ) > 0)
    stop ("sample names do not match the colnames of count table ")
  if ( sum ( sample_annot$sample_names != rownames(design) ) > 0)
    stop ("sample names do not match the rownames of desing matrix ")
  if ( sum( colnames(design) != colnames(contrasts) ) > 0)
    stop ("colnames(design) do not match the colnames of contrasts")
  if (sum( rownames(design) != colnames(counts) ) > 0)
    stop ("rownames(design) do not match the colnames of counts")

  dge <- edgeR::DGEList(counts,samples = sample_annot)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design = design, robust=robust)
  fit <- edgeR::glmQLFit(dge, design = design , robust=robust)

  if(plots==TRUE){
    edgeR::plotBCV(dge)
    edgeR::plotQLDisp(fit)
  }
  de_results = list()
  for(i in 1:nrow(contrasts)){
    de_table = edgeR::glmQLFTest(glmfit = fit, contrast = as.numeric(  contrasts[i,] ) )$table
    de_table$F = NULL
    colnames(de_table) = c( "log2FoldChange", "logCPM", "pvalue")
    de_table$padj = p.adjust(de_table$pvalue , method = "BH")
    de_table$DE = 0
    de_table$DE [de_table$padj < adj_pv_threshold &
                   de_table$log2FoldChange >  log2fc_threshold] =  1
    de_table$DE [de_table$padj < adj_pv_threshold &
                   de_table$log2FoldChange < -log2fc_threshold] = -1

    de_results[[i]] = de_table
  }
  names(de_results) = rownames(contrasts)

  return(de_results)
}


#' DeSeq2 wrapper
#' @param counts data.matrix with the counts data, rownames should represent features and colnames sample names
#' @param sample_annot data.frame with sample annotation sample_annot$sample_names should be same as colnames(counts), it has to contain exp_fac column with experimental factor, control needs to be a the first level of the factor (see relevel() )
#' @param adj_pv_threshold adjusted p-value significance threshold for DE testing, 0.05 by default
#' @param log2fc_threshold log2 Fold Change threshold for DE testing, 1 by default
#' @return a list of differential expression results, each element of the list corresponds to one comparison: each experimental factor level vs control (first factor level)
#' @importFrom  DESeq2 DESeqDataSetFromMatrix DESeq resultsNames lfcShrink
#' @export
deseq_DE_wrapper = function( counts, sample_annot, adj_pv_threshold = 0.05, log2fc_threshold = 1){

  #initial checks
  if ( sum ( sample_annot$sample_names != colnames(counts) ) > 0)
    stop ("sample names do not match the colnames of count table ")
  if (is.null(sample_annot$exp_fac))
    stop ("exp_fac column not found in sample_annot")

  #construction DESeq object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = sample_annot,
                                design    = ~ exp_fac)
  # run DESeq
  dds <- DESeq(dds)

  comparisons = resultsNames(dds)
  comparisons = comparisons [ -grep ("Intercept" , comparisons)]

  # Shrink log2 fold changes using apeglm
  de_results = list()
  for( comparison in comparisons) {
    de_results [[ comparison ]] <- lfcShrink(dds = dds, coef=comparison, type="apeglm")
  }

  for (i in 1:length(de_results)){
    de_results[[i]]$padj [ is.na(de_results[[i]]$padj)] = 1
    de_results[[i]]$DE = 0
    de_results[[i]]$DE [de_results[[i]]$padj < adj_pv_threshold &
                          de_results[[i]]$log2FoldChange >  log2fc_threshold] =  1
    de_results[[i]]$DE [de_results[[i]]$padj < adj_pv_threshold &
                          de_results[[i]]$log2FoldChange < -log2fc_threshold] = -1
  }
  names(de_results) = sub("exp_fac_" , "" ,  names(de_results))

  de_results
}



#' limma differential analysis wrapper
#' @param data a matrix of expression values or activities, data should be normally distributed (e.g. normalized and log2 transformed)
#' @param design design matrix
#' @param contrasts matrix with contrasts
#' @param adj_pv_threshold adjusted p-value significance threshold for DE testing, 0.05 by default
#' @param log2fc_threshold log2 Fold Change threshold for DE testing, 1 by default#'
#' @return a list of differential expression results, each element of the list corresponds to one row of contrast matrix
#' @export
#' @importFrom  limma lmFit contrasts.fit eBayes topTable
limma_DE_wrapper = function( data , design, contrasts, log2fc_threshold = 1,adj_pv_threshold = 0.05 ){
  fit  = limma::lmFit( data , design )
  fit2 = limma::contrasts.fit(fit, t(contrasts))
  fit2 = limma::eBayes(fit2)

  de_results = list()
  for(i in 1:nrow(contrasts)){
    de_table = limma::topTable(fit2, coef=i, adjust="BH", number = nrow(data), sort.by = "none")
    de_table$B = NULL
    de_table$t = NULL
    colnames(de_table) = c( "log2FoldChange", "Average", "pvalue", "padj")
    de_table$padj = p.adjust(de_table$pvalue , method = "BH")
    de_table$DE = 0
    de_table$DE [de_table$padj < adj_pv_threshold &
                   de_table$log2FoldChange >  log2fc_threshold] =  1
    de_table$DE [de_table$padj < adj_pv_threshold &
                   de_table$log2FoldChange < -log2fc_threshold] = -1

    de_results[[i]] = de_table
  }
  names(de_results) = rownames(contrasts)
  return(de_results)
}


