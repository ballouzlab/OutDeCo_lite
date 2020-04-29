#' Reformat DE output
#'
#' @param output DE output
#' @param method method used
#' @return \code{combined_out}
#' @examples
#'
#' @import edgeR DESeq2 stats
#' @export
#'
#'

reformat_degs <- function(output, method = "wilcox") {

    if (method == "edgeR") {
        padj <- p.adjust(output$table$PValue )
        degs <- data.frame(  mean_cpm = 2^output$table$logCPM ,
                             log2_fc = output$table$logFC,
                             pvals = output$table$PValue,
                             padj = padj )
    }
    if (method == "DESeq2") {
        degs <- data.frame(  mean_cpm = output$baseMean,
                             log2_fc = output$log2FoldChange,
                             pvals = output$pvalue,
                             padj =  output$padj )

    }
    combined_out <- list(degs, output)
    names(combined_out) = c("degs", "output")
    return( combined_out )
}
