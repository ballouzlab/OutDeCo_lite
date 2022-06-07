#' Gene set enrichment
#'
#' @param genes an array of numeric values
#' @param genes.labels GO matrix or some such
#' @param voc term descriptions
#' @return \code{results}
#' @examples
#' data <-  rnorm(1000)
#' geo_mean(data)
#'
#' @import stats utils
#' @export


gene_set_enrichment <- function(genes, genes.labels, voc){

  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p

  m = match ( genes, genes.names )
  filt.genes  = !is.na(m)
  filt.labels = m[filt.genes]


  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g

  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]

  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
  else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets

  test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  sigs = pvals < ( 0.05/length(pvals) )
  pvals.adj = p.adjust( pvals, method="BH")

  results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj", "sig" )

  return (results)

}




#' Gene set enrichment - AUC version
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


