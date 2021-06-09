
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
