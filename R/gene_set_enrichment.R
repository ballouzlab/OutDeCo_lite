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
