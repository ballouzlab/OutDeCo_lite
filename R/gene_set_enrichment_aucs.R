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

o = order(gene_rankings)
nbins = round(length(o)/100)

enrich_mat = sapply(1:length(gene_set_aucs), function(j) (sapply(0:nbins, function(i) sum(gene_sets[o,j][ (1:100)+ i*100 ]  )  )  ) )  
enrich_mat[is.na(enrich_mat)] = 0 
enrich_mat = t(enrich_mat)/colSums(enrich_mat)
#enrich_mat = sapply(1:length(gene_set_aucs), function(j) (sapply(0:nbins, function(i) sum(gene_sets[o,j][ (1:1000)+ i*1000 ]  )  )  ) )  
heatmap.2(  (enrich_mat), col=colssig, trace="none", density="none", Rowv=F, Colv=F)

  return (gene_set_aucs)

}


