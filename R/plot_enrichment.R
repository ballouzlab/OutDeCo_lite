#' Plot recurrence histogram
#'
#' @param data values
#' @param n_max total number of studies
#' @examples
#'
#' data <-  matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' n_max <- dim(data)[2]
#' plot_enrichment( data, n_max)
#'
#' @import viridis venn
#' @importFrom gplots colorpanel
#' @export

plot_enrichment <- function(data, n_max, ...) {

   labels <- names(data)
   paths <-  sapply(1:n_max, function(i) (data[[i]]$padj<0.05)*1 )
   paths_padj <-  sapply(1:n_max, function(i) (data[[i]]$padj) )
   rownames(paths) <-  data[[1]][,1]
   rownames(paths_padj) <- data[[1]][,1]
   colnames(paths) <- labels
   colnames(paths_padj) <- labels


   f <- rowSums( paths) > 1
   sigtemp <- paths
   sigtemp[sigtemp==1] <- "*"
   sigtemp[sigtemp==0] <- ""
   colssig <- colorpanel(100, "white", "red", "darkmagenta")

   heatmap.2( -log10(t(paths_padj[f,])), Rowv=F,
              col=colssig, cexRow = 0.7, cexCol = 0.7,
              cellnote=t(sigtemp[f,]),
              notecol="black",
              notecex=1,
              keysize=1,
              key.xlab="-log10 adjusted P-value",
              key.title="Enrichment", trace="none", density="none" )

}




