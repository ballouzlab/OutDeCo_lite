#' Plot recurrence histogram
#'
#' @param data enrichment results
#' @param gene_list genes
#' @param gene_sets gene sets
#' @import viridis venn
#' @importFrom gplots colorpanel
#' @export

plot_gene_set_enrichment <- function(data, gene_list, gene_sets) {
#cols4 = colorpanel(63, "lightgrey", "blue", "darkblue")
#cols5 = colorpanel(300, "lightgrey", "red", "darkred")
colssig = colorpanel(100, "lightgrey", "red", "darkmagenta")

   paths <- data$term
   paths_padj <- data$padj
   paths_pval <- data$pval

   # Match gene list with genes in gene sets 
   m <-  match( gene_list, rownames(gene_sets) ) 
   sub_gene_sets <- gene_sets[m,]
   sub_gene_sets[is.na(sub_gene_sets)] = 0
   rownames(sub_gene_sets) <- gene_list

   # Match enrichment results with gene sets input 
   m <- match( paths, colnames(gene_sets) )
   sub_gene_sets <- sub_gene_sets[,m]

   # Set up pval and padj matrices to plot
   sub_gene_sets_pval <- (sub_gene_sets * 0 ) + 1 
   sub_gene_sets_padj <- (sub_gene_sets * 0 ) + 1 

   for(i in 1:length(paths_pval)){ 
       sub_gene_sets_pval[ sub_gene_sets[,i]==1,i] = paths_pval[i]
       sub_gene_sets_padj[ sub_gene_sets[,i]==1,i] = paths_padj[i]
   }

   log10sub_gene_sets_pval = -log10(sub_gene_sets_pval)
   log10sub_gene_sets_padj = -log10(sub_gene_sets_padj)
   log10sub_gene_sets_pval[!is.finite(log10sub_gene_sets_pval)] = 0 
   log10sub_gene_sets_padj[!is.finite(log10sub_gene_sets_padj)] = 0 
   
    filt = colSums(log10sub_gene_sets_pval) > 0

   #sigtemp <- 1*(sub_gene_sets_padj < 0.05)
   #sigtemp[sigtemp==1] <- "*"
   #sigtemp[sigtemp==0] <- ""

   heatmap.2(log10sub_gene_sets_pval[,filt]  , Rowv=F, Colv=F,
              col=colssig, cexRow = 0.7, cexCol = 0.7,
    #          cellnote=sigtemp,
              notecol="black",
              notecex=1,
              keysize=1,
              key.xlab="-log10 adjusted P-value",
              key.title="Enrichment", trace="none", density="none" )
      
}




