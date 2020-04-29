#' Probability calc
#'
#' @param gene_sets A matrix of genes (rows)  by study (columns)
#' @return \code{fdrs_bin}
#' @examples
#'
#' gene_sets <- matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' fdrs_bin <- calc_binom_recur(gene_sets)
#'
#' @import stats utils
#' @export
#'

calc_binom_recur <- function(gene_sets) {
    n_sets <- dim(gene_sets)[2]
    n_genes <- dim(gene_sets)[1]
    p <- max(colSums(gene_sets)) / n_genes
    fdrs_bin <- p.adjust(pbinom(0:n_sets, n_genes, p, lower.tail = F), n = n_genes)
    fdrs_bin <- cbind(fdrs_bin, 0:n_sets)
    return(fdrs_bin)
}
