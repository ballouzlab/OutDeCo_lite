#' Gene set enrichment
#'
#' @param gene_sets GO matrix or some such
#' @param gene_rankings an array of numeric values
#' @return \code{results}
#' @importFrom EGAD auc_multifunc
#' @export


gene_set_enrichment_aucs <- function(gene_sets, gene_rankings){

  
  m <- match( rownames(gene_sets), names(gene_rankings) )
  f.g = !is.na(m)
  f.r = m[f.g]
  gene_sets = gene_sets[f.g,]
  gene_rankings = rank(gene_rankings[f.r])

  gene_set_aucs <- EGAD::auc_multifunc(gene_sets, gene_rankings )
   names(gene_set_aucs) = colnames(gene_sets)
  return (gene_set_aucs)

}
