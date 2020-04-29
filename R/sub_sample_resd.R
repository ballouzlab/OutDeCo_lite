#' Connectivity test wrapper function
#'
#' @param net A matrix (network)
#' @param dir output directory
#' @param label name of study
#' @param gene_list Gene set to test
#' @param n_sub minimum gene set size of studies being compared
#' @param n_r number of repeats
#' @param filename fileanme
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)
#' studies <- paste0('study', 1:10 )
#' studies_genes   <-  paste0('gene', 1:100 )
#' colnames(genesets) <- studies
#' rownames(genesets) <- studies_genes
#' i <- 3
#' gene_list <- as.matrix(cbind(studies_genes,
#'              as.numeric(genesets[,i]))[genesets[, i] > 0, ])[,1]
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:1000 )
#' colnames(network) <-paste0('gene', 1:1000 )
#'
#' nettype <- 'random'
#' dir <- 'out_test'
#' label <- studies[i]
#' filename <- 'output'
#'
#' # sub_sample_resd(network, dir, label, gene_list, n_sub, n_r, filename)
#'
#' @import stats utils
#' @export
#'
#'

sub_sample_resd <- function(net, dir, label, gene_list, n_sub, n_r, filename) {

    if (dim(gene_list)[1] < n_sub) {
        residual_connectivity_score(net, dir, label, gene_list)
    } else {
        residual_connectivity_score_sub(net, dir,
                                        label, gene_list,
                                        n_sub, n_r, filename)

    }
}
