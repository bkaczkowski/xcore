#' Fisher's test signatures
#'
#' Perform Fisher test between binary matrices columns. The function can be work
#' using parallelization  if  a parallel backend is registered (eg. \code{doMC}).
#'
#' @inheritParams stats::fisher.test
#' @param exprss_mat a binary expression matrix.
#' @param sign_mat a binary signature matrix.
#'
#' @return a list with \code{"p.value"} matrix of p-values and \code{"odds"}
#'   matrix of odds ratios from the Fisher test.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom stats stats
#'
#' @export
fisherTestSignature <- function(exprss_mat, sign_mat, alternative = "greater") {
  stopifnot(is.matrix(exprss_mat)) # || is(exprss_mat, "dgCMatrix")) #TODO if we support sparseMatrices all methods eg. colSums must also support it
  stopifnot(is.matrix(sign_mat)) # || is(sign_mat, "dgCMatrix"))
  stopifnot(nrow(exprss_mat) == nrow(sign_mat))
  
  ecs <- colSums(exprss_mat)
  if (any(drop <- ecs == nrow(exprss_mat) | ecs == 0)) {
    warning(sprintf("Droping %i constant columns from expression", sum(drop)))
    exprss_mat <- exprss_mat[, ! drop]
  }
  esi <- colSums(sign_mat)
  if (any(drop <- esi == nrow(sign_mat) | esi == 0)) {
    warning(sprintf("Droping %i constant columns from signatures", sum(drop)))
    sign_mat <- sign_mat[, ! drop]
  }
  
  fishers <- foreach(j = seq_len(ncol(exprss_mat))) %dopar%
    apply(X = sign_mat,
          MARGIN = 2,
          FUN = function(sig) fisher.test(x = exprss_mat[, j, drop = TRUE],
                                          y = sig,
                                          alternative = alternative)
    )
  
  .getVal <- function(fishers, val, colnm) {
    mat <- do.call(what = cbind,
                   args = lapply(X = fishers,
                                 FUN = function(x) {
                                   do.call(what = rbind,
                                           args = lapply(x, `[[`, val))
                                 }))
    colnames(mat) <- colnm
    mat
  }
  
  list(
    p.value = .getVal(fishers, "p.value", colnames(exprss_mat)),
    odds = .getVal(fishers, "estimate", colnames(exprss_mat))
  )
}

#' Normalize the count table with edgeR
#'
#' @param counts expression table with counts
#' @param method passed to to edgeR::calcNormFactors see ?edgeR::calcNormFactors
#' @param log passed on to edgeR::cpm , should the normalized counts be log2 transformed (default FALSE)
#' @param prior.count average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
#' @param ... other parameter passed on to edgeR::calcNormFactors or edgeR::cpm
#' @return normalized count table
#' @import edgeR
#' @export
normalize_counts = function(counts,
                            method = "RLE",
                            log = FALSE ,
                            prior.count = 1 ,
                            ...) {
  dge = edgeR::DGEList(counts)
  dge = edgeR::calcNormFactors(dge , method = method , ...)
  cpm = edgeR::cpm(dge , log = log, prior.count = prior.count, ...)
  cpm
}

#' Select top variable rows
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_var = function(matrix, n = 30) {
  act_var <- apply(matrix, 1, sd)
  matrix <- matrix[order(act_var, decreasing = TRUE), ]
  matrix[1:n, ]
}

#' Select rows with highest mean
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_mean_up = function(matrix,  n = 30) {
  act_mean <- apply(matrix, 1, mean)
  matrix <- matrix[order(act_mean, decreasing = TRUE), ]
  matrix[1:n, ]
}

#' Select rows with lowest mean
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_mean_down = function(matrix,  n = 30) {
  act_mean <- apply(matrix, 1, mean)
  matrix <- matrix[order(act_mean, decreasing = FALSE), ]
  matrix[1:n, ]
}

#' Annotation Worker by Same Strand Overlap
#'
#' Usually used by the gencode_one_direction_annotator function.
#' @param regions query genomic regions as GRanges object
#' @param gencode Gencode annotation as GRanges object, usually of one type like gencode$type == "exon"
#' @param annotation_label label (character) to mark what annotation was used, e.g. "gene", "exon" or "promoters"
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges findOverlaps
#' @export
gencode_one_direction_findOverlaps = function( regions , gencode, annotation_label = " "){
  # ordering so higher level gene annotation will overwrie lower level annotation
  # if multiple genes are annotated to same region
  gencode = gencode[ order( gencode$level , decreasing = T),]
  
  if( is.null(regions$annotation) )   { regions$annotation = ""}
  if( is.null(regions$symbol) )       { regions$symbol       = ""}
  if( is.null(regions$gene_type) )    { regions$gene_type    = ""}
  if( is.null(regions$level) )        { regions$level        = ""}
  
  hits = GenomicRanges::findOverlaps( regions, gencode , ignore.strand=FALSE)
  regions$symbol[hits@from]       = gencode$gene_name[hits@to]
  regions$gene_type[hits@from]    = gencode$gene_type[hits@to]
  regions$level[hits@from]        = gencode$level[hits@to]
  regions$annotation[hits@from]   = annotation_label
  
  regions
}

#' Annotating Same Strand by Overlap
#' A wrapper that calls gencode_one_direction_findOverlaps on all types of Gencode Annotation
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @export
gencode_one_direction_annotator = function( regions , gencode ){
  
  gencode_gene        = gencode[ gencode$type == "gene"]
  gencode_exon        = gencode[ gencode$type == "exon"]
  gencode_3_prime_UTR = gencode[ gencode$type == "three_prime_UTR"]
  gencode_5_prime_UTR = gencode[ gencode$type == "five_prime_UTR"]
  gencode_promoters   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 500, downstream = 500, use.names=TRUE)
  
  regions = gencode_one_direction_findOverlaps( regions , gencode_gene        , annotation_label = "intron" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_exon        , annotation_label = "exon" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_3_prime_UTR , annotation_label = "three_prime_UTR" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_5_prime_UTR , annotation_label = "five_prime_UTR" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_promoters   , annotation_label = "promoter" )
  
  regions
}

#' Annotating Both Strands by Overlap
#'
#' A wrapper that calls gencode_one_direction_annotator on same and inverted strand
#' to annotate overlap on both strands
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges strand invertStrand
#' @export
gencode_two_direction_annotator = function(regions , gencode  ){
  
  regions_inverted_strand = regions
  GenomicRanges::strand(regions_inverted_strand) = GenomicRanges::invertStrand(GenomicRanges::strand(regions))
  
  regions  = gencode_one_direction_annotator( regions = regions, gencode = gencode)
  regions_inverted_strand  = gencode_one_direction_annotator( regions = regions_inverted_strand, gencode = gencode)
  
  regions$opposite_strand_annotation = regions_inverted_strand$annotation
  regions$opposite_strand_symbol = regions_inverted_strand$symbol
  regions$opposite_gene_type = regions_inverted_strand$gene_type
  
  regions
}

#' Annotating Strandless regions by Overlap
#'
#' This is mostly for bidirectional/strandless enhancer regions
#' The wrappers call the gencode_one_direction_annotator twice after setting the strand to + and -.
#' Finally the output is nicely summarized
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges strand mcols
#' @export
gencode_strandless_annotator = function( regions, gencode , output_option = 2) {
  
  regions_tmp = regions
  GenomicRanges::strand(regions_tmp) = "+"
  GenomicRanges::mcols(regions_tmp)= NULL
  regions_annotated_plus  = gencode_one_direction_annotator( regions = regions_tmp, gencode = gencode)
  regions_annotated_plus$level = NULL
  regions_annotated_plus  = GenomicRanges::mcols(regions_annotated_plus)
  
  GenomicRanges::strand(regions_tmp) = "-"
  regions_annotated_minus = gencode_one_direction_annotator( regions = regions_tmp, gencode = gencode)
  regions_annotated_minus$level = NULL
  regions_annotated_minus = GenomicRanges::mcols(regions_annotated_minus)
  
  if( output_option == 1 ){
    regions$strandless_annotation = paste( regions_annotated_plus$annotation, "(+);", regions_annotated_minus$annotation , "(-)", sep = "" )
    regions$strandless_symbol = paste( regions_annotated_plus$symbol, "(+);", regions_annotated_minus$symbol , "(-)", sep = "" )
    regions$strandless_gene_type = paste( regions_annotated_plus$gene_type, "(+);", regions_annotated_minus$gene_type , "(-)", sep = "" )
    
  }else if( output_option == 2) {
    colnames(regions_annotated_minus) = paste( colnames(regions_annotated_minus), "_minus",sep = "")
    colnames(regions_annotated_plus)  = paste( colnames(regions_annotated_plus ), "_plus" ,sep = "")
    GenomicRanges::mcols(regions) = DataFrame( cbind(GenomicRanges::mcols(regions) ,  regions_annotated_plus, regions_annotated_minus))
  }
  
  regions
}

#' Annotation Worker by Same Strand DistanceToNearest
#'
#' @param regions query genomic regions as GRanges object
#' @param gencode Gencode annotation as GRanges object, usually of one type like gencode[ gencode$type == "exon"]
#' @param annotation_label label (character) to mark what annotation was used, e.g. "gene", "exon" or "promoters"
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges distanceToNearest
#' @export
gencode_one_direction_distanceToNearest = function( regions , gencode){
  
  # ordering so higher level gene annotation will overwrie lower level annotation
  # if multiple genes are annotated to same region
  gencode = gencode[ order( gencode$level , decreasing = F),]
  
  if( is.null(regions$nearest_symbol) )       { regions$nearest_symbol         = ""}
  if( is.null(regions$nearest_gene_type) )    { regions$nearest_gene_type      = ""}
  if( is.null(regions$nearest_distance) )     { regions$nearest_distance       = -1}
  
  hits = GenomicRanges::distanceToNearest( regions, gencode , ignore.strand=FALSE)
  regions$nearest_symbol[hits@from]       = gencode$gene_name[hits@to]
  regions$nearest_gene_type[hits@from]    = gencode$gene_type[hits@to]
  regions$nearest_distance                = hits@elementMetadata$distance
  
  regions
}

#' Annotate with Nearest Promoter for Both Strands
#'
#' The wrappers call the gencode_one_direction_distanceToNearest twice after setting the strand to + and -.
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "gencode.v31.annotation.gff3.gz") ,  (gencode[ gencode$type == "transcript"] is also OK)
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges promoters mcols strand
#' @importFrom S4Vectors DataFrame
#' @export
gencode_nearest_promoter_both_strands = function( regions, gencode ) {
  
  gencode_promoters   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 0, downstream = 0, use.names=TRUE)
  
  regions_tmp = regions
  GenomicRanges::mcols(regions_tmp)= NULL
  
  GenomicRanges::strand(regions_tmp) = "+"
  nearest_promoter_plus  = gencode_one_direction_distanceToNearest( regions = regions_tmp, gencode = gencode_promoters)
  nearest_promoter_plus  = GenomicRanges::mcols(nearest_promoter_plus)
  colnames(nearest_promoter_plus) = paste( colnames(nearest_promoter_plus), "_plus",sep = "")
  
  GenomicRanges::strand(regions_tmp) = "-"
  nearest_promoter_minus = gencode_one_direction_distanceToNearest( regions = regions_tmp, gencode = gencode_promoters)
  nearest_promoter_minus = GenomicRanges::mcols(nearest_promoter_minus)
  colnames(nearest_promoter_minus) = paste( colnames(nearest_promoter_minus), "_minus",sep = "")
  
  GenomicRanges::mcols(regions) = DataFrame( cbind( GenomicRanges::mcols(regions) , nearest_promoter_plus, nearest_promoter_minus))
  
  regions
}

#' Annotate with Nearest Promoter on the Same Strand
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "gencode.v31.annotation.gff3.gz") ,  (gencode[ gencode$type == "transcript"] is also OK)
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges promoters mcols
#' @importFrom S4Vectors DataFrame
#' @export
gencode_nearest_promoter_same_strand = function( regions, gencode , cut_off_distance = 500) {
  
  gencode_promoters   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 0, downstream = 0, use.names=TRUE)
  
  regions = gencode_one_direction_distanceToNearest( regions = regions, gencode = gencode_promoters)
  
  if( !is.null( cut_off_distance )){
    too_far = regions$nearest_distance > cut_off_distance
    regions$nearest_symbol     [ too_far ]  = ""
    regions$nearest_gene_type  [ too_far ]  = ""
  }
  
  regions
}

#' Simplify Gencode Annotation
#' @param annotation character vector with Gencode annotation that can be used to grep() labels as: "intron", "exon", "three_prime_UTR", "five_prime_UTR", "promoter"
#' @return character vector with simplified labels
#' @export
simplify_gencode_annotation = function( annotation ) {
  simplified_annotation = rep ( "intergenic", length(annotation))
  simplified_annotation [ grep( "intron" , annotation)] = "intron"
  simplified_annotation [ grep( "exon"   , annotation)] = "exon"
  simplified_annotation [ grep( "three_prime_UTR" , annotation)] = "three_prime_UTR"
  simplified_annotation [ grep( "five_prime_UTR"  , annotation)] = "five_prime_UTR"
  simplified_annotation [ grep( "promoter"        , annotation)] = "promoter"
  simplified_annotation
}

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

#' Read Collections of MSigDB signatures
#' @param path_to_gmts a path to directory containing .gmt files downloaded from http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#' @param msigdb_version a version of MSigDB to read in, e.g. ".v7.0.entrez.gmt"
#' @return a list, each element of the list represents one signature, e.g. c2.cgp___RUNNE_GENDER_EFFECT_UP" represent "RUNNE_GENDER_EFFECT_UP" signature from "c2.cgp" collection
#' @importFrom fgsea gmtPathways
#' @export
read_msigdb_gmts = function( path_to_gmts , msigdb_version = ".v7.0.entrez.gmt" ){
  msigdb = data.frame ( file = list.files( path_to_gmts, pattern = msigdb_version, full.names = T) , stringsAsFactors = F)
  msigdb$collection = gsub(".*/", "" , msigdb$file)
  msigdb$collection = gsub(msigdb_version, "" , msigdb$collection)
  gene_sets = list()
  for (i in 1:nrow(msigdb) ){
    gene_set <- fgsea::gmtPathways(msigdb$file[i])
    names( gene_set ) = paste( msigdb$collection[i] , names( gene_set ), sep = "___")
    gene_sets[[i]] = gene_set
  }
  gene_sets = do.call(gene_sets, what = c)
  gene_sets = gene_sets [! duplicated ( sub(".*___" , "", names(gene_sets) ) ) ]
  gene_sets
}

#' Rank Genes for fgsea
#'
#' @param de_results data.frame with DE results from edgeR or DeSeq, should contail "pvalue" and "log2FoldChange" columns and DE stats for all features
#' @param based_on choose if "pvalue" or "log2FoldChange" should be used to rank the genes/featrues
#' @return sorted vector of ranks with names that correspond to gene names (required input to fgse() function)
#' @export
rank_genes = function( de_results , based_on = c("pvalue" , "log2FoldChange")){
  if( based_on == "pvalue"){
    pv_ranks = -log10(de_results$pvalue) * sign(de_results$log2FoldChange)
    names(pv_ranks) = rownames(de_results)
    return(sort(pv_ranks))
  }else if ( based_on == "log2FoldChange")  {
    fc_ranks = de_results$log2FoldChange
    names(fc_ranks) = rownames(de_results)
    return( sort(fc_ranks) )
  } else {
    stop ("based_on parameter not recognized")
  }
}

#' Get Matrix of Fold Changes (FCs)
#'
#' @param de_res a list of DE results data frame with log2FoldChange column
#' @return matrix of fold changes
#' @export
get_fc_matrix = function( de_res){
  extract_fc = function( x ) { x$log2FoldChange }
  fc_list = lapply( de_res, extract_fc)
  fc_matrix = do.call(what = cbind, fc_list)
  rownames(fc_matrix) = rownames(de_res[[1]])
  return(fc_matrix)
}

#' Extract a column of interest from list of data.frame's
#' @param de_res list of data.frame objects, e.g. DE results
#' @param column name of the column to extract, e.g. "log2FoldChange"
#' @return extracted column cbind'ed into a matrix
#' @export
extract_a_column = function( de_res , column ){
  extract_fun = function( x ) { x [ , column ] }
  col_list = lapply( de_res, extract_fun)
  col_matrix = do.call(what = cbind, col_list)
  rownames(col_matrix) = rownames(de_res[[1]])
  return(col_matrix)
}

#' Create a list of signatures from 0,1 encoded matrix
#' @param int_mat 0,1 encoded matrix with columns representing signatures and rows representing the genes, 1 means that a gene belongs to a pathway
#' @return a list gene signatures
#' @export
matrix_to_list = function( int_mat ){
  sig_list = list()
  for ( j in 1:ncol(int_mat) ) {
    signature_elements = int_mat [,j] == 1
    sig_list[[j]] = rownames(int_mat) [ signature_elements ]
  }
  names(sig_list) = colnames(int_mat)
  sig_list
}

#' Create a 0,1 encoded matrix from a list of signatures
#' @param sig_list a list gene signatures
#' @return 0,1 encoded matrix with columns representing signatures and rows representing the genes, 1 means that a gene belongs to a pathway
#' @export
list_to_matrix = function( sig_list ) {
  unique_features = sort( unique ( do.call(c, msigdb) ) )
  overlap_mat = matrix (0L, nrow = length(unique_features), ncol = length(msigdb) )
  colnames(overlap_mat) = names(msigdb)
  rownames(overlap_mat) = unique_features
  for ( j in 1:ncol(overlap_mat) ) {
    overlap_mat[ unique_features  %in% msigdb[[j]]  , j ] = 1L
  }
  overlap_mat
}

#' Create a sparseMatrix from a list of signatures
#' @param sig_list a list gene signatures
#' @return a sparseMatrix, columns represent signatures and rows represent the genes, TRUE means that a gene belongs to a pathway
#' @export
#' @importFrom Matrix sparseMatrix
list_to_sparse_matrix = function( sig_list ) {
  unique_features = sort( unique ( do.call(c, msigdb) ) )
  overlap_mat = Matrix::sparseMatrix (dims = c(length(unique_features), length(msigdb) ), i={}, j={} )
  colnames(overlap_mat) = names(msigdb)
  rownames(overlap_mat) = unique_features
  for ( j in 1:ncol(overlap_mat) ) {
    overlap_mat[ unique_features  %in% msigdb[[j]]  , j ] = TRUE
  }
  overlap_mat
}


#' Create list of ranks from table of FoldChanges
#' @param fc_table a table of FoldChanges
#' @return list of rank vectors, each rank vector can be used as input to fgsea() function
#' @export
ranks_from_FC_table = function(fc_table){
  ranks = list()
  for ( j in 1:ncol( fc_table ) ) {
    rank = fc_table[, j]
    names(rank) = rownames(fc_table)
    ranks[[j]] = sort(rank)
  }
  names(ranks) = colnames(fc_table)
  ranks
}

#' Run fgsea() wrapper/convenience function
#'
#' @param ranks a list of ranks as produced by ranks_from_FC_table()
#' @param pathways a list of signatures as produced by matrix_to_list()
#' @param nperm number of permutation done by fgsea(), default nperm = 1 to calculate the enriment scores without calculating p-value, you can calculate p-value by using nperm = 1000 or nperm = 10000 but it can be prohibitively slow
#' @return a list, each element represent fgsea() result for one differential expression comparison, across all signatures
#' @importFrom fgsea fgsea
#' @export
run_fgsea = function( ranks , pathways, nperm = 1){
  fgseaRes = list()
  for ( rank in   names(ranks) ) {
    fgseaRes [[rank]] <- fgsea::fgsea(pathways = pathways, stats = ranks[[rank]],
                                      minSize=1, maxSize=50000, nperm=nperm)
  }
  fgseaRes
}
#' Make Enrichment Score Matrix
#'
#' @param fgseaRes list of fgsea() results as generated by run_fgsea()
#' @param score choose between two enrichment scores: "NES" or "ES"
#' @return a matrix with enrichment scores ("NES" or "ES"), each column represents DE comparison, each row represents a pathway/signature
#' @export
make_enrichment_score_matrix= function( fgseaRes, score = c("NES", "ES") ) {
  FUN = function(x) {x[[score]]}
  scores = do.call( cbind , ( lapply(  fgseaRes, FUN )))
  rownames( scores ) = fgseaRes[[1]]$pathway
  scores [ is.nan(scores) ]  = 0
  scores
}

#' Calculate Activity Scores (as Mean FC of Singature Features)
#' @param signature_mat a sparseMatrix, columns represent signatures and rows represent the genes, 1 means that a gene belongs to a pathway
#' @param de_res list of data.frame objects, e.g. DE results
#' @param min_genes_per_sig minimum number of genes per signature, signatures with fewer genes are discarded
#' @return a matrix of activities, rows represent signatures and columns represent the DE comparisons
#' @export
#' @importFrom Matrix t colSums crossprod
de_res_to_activity_scores = function( signature_mat ,de_res , min_genes_per_sig = 10 ) {
  fc_matrix = xcore::get_fc_matrix( de_res )
  fc_matrix [ is.na (fc_matrix)] = 0
  fc_matrix = fc_matrix [ rownames(fc_matrix) %in% rownames( signature_mat ) , ]
  fc_matrix = fc_matrix [ order( rownames(fc_matrix) ) , ]
  fc_matrix = scale(fc_matrix)
  fc_matrix = data.matrix( fc_matrix  )
  
  signature_mat = signature_mat [ rownames(signature_mat) %in% rownames(fc_matrix), ]
  signature_mat = signature_mat [ order( rownames(signature_mat) ) , ]
  signature_mat = signature_mat [ , Matrix::colSums(signature_mat) >= min_genes_per_sig]
  
  signature_mat_scaled  = Matrix::t(Matrix::t(signature_mat) / Matrix::colSums(signature_mat ))
  
  if ( sum (rownames(fc_matrix) != rownames(signature_mat_scaled) ) > 0 ) {
    stop ("DE results and count table features don't match")
  }
  Matrix::crossprod( signature_mat_scaled ,fc_matrix  )
}

#' Calculate Activity Scores for Expressim Matrix
#' @param signature_mat sparseMatrix, columns represent signatures and rows represent the genes, 1 means that a gene belongs to a pathway
#' @param counts a matrix with expression counts
#' @param min_genes_per_sig minimum number of genes per signature, signatures with fewer genes are discarded
#' @param prior.count number of counts to be added as an offset in log2 transformation, used by edgeR::cpm()
#' @param control_samples vector of indices of control samples to be used as reference, NULL by default
#' @return a matrix of activities, rows represent signatures and columns represent the samples
#' @export
#' @importFrom Matrix t colSums crossprod
counts_to_activity_scores = function( signature_mat ,counts , min_genes_per_sig = 10 , prior.count = 4 , control_samples = NULL) {
  
  signature_mat = signature_mat [ rownames(signature_mat) %in% rownames(counts), ]
  signature_mat = signature_mat [ order( rownames(signature_mat) ) , ]
  signature_mat = signature_mat [ , Matrix::colSums(signature_mat > 0) >= min_genes_per_sig]
  signature_mat_scaled  = Matrix::t(Matrix::t(signature_mat) / Matrix::colSums(signature_mat ))
  
  tpm = xcore::normalize_counts(counts, method = "RLE", prior.count = prior.count, log = TRUE)
  tpm = tpm [ rownames(tpm) %in% rownames(signature_mat_scaled) , ]
  tpm = tpm [ order( rownames(tpm) ) , ]
  tpm = data.matrix( tpm  )
  
  if( is.null (control_samples) ){
    tpm = tpm - rowMeans(tpm)
  }else{
    tpm = tpm - rowMeans(tpm[ ,  control_samples ])
  }
  
  if ( sum (rownames(tpm) != rownames(signature_mat_scaled) ) > 0 ) {
    stop ("Signature mat and count table features don't match")
  }
  Matrix::crossprod( signature_mat_scaled , tpm  )
}
