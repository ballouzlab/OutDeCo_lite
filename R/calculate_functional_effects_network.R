#' Connectivity tests
#'
#' @param genesets Sets of genes to test
#' @param network A matrix (network)
#' @param nettype network label
#' @param dir output directory
#' @param studies gene set labels
#' @param nsub minimum gene set size of studies being compared (downsample to this)
#' @param nr number of repeats
#' @return \code{output} list with effect sizes for network connectivities
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)
#' studies <- paste0('study', 1:10 )
#' genes   <-  paste0('gene', 1:100 )
#' colnames(genesets) <- studies
#' rownames(genesets) <- genes
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:1000 )
#' colnames(network) <-paste0('gene', 1:1000 )
#' nettype <- 'random'
#' dir <- 'out_test'
#'
#'
#' @import stats utils
#' @export
#'
# calculate_functional_effects_network(genesets, network, nettype, dir, studies, 5)

calculate_functional_effects_network <- function(genesets, network, nettype, dir, studies, nsub, nr = 1000) {

    studies_genes <- rownames(genesets)
    gene_lists <- lapply(1:length(studies), function(i) as.matrix(cbind(studies_genes, as.numeric(genesets[, i]))[genesets[,
        i] > 0, ]))
    labels <- lapply(1:length(studies), function(i) paste(studies[i], nettype, sep = "."))
    file <- paste(dir, "/", "random.", nsub, ".", nettype, ".residuals.random", sep = "")

    an <- lapply(1:length(studies), function(i) residual_connectivity_score(network, dir, labels[[i]], gene_lists[[i]]))
    an <- matrix(unlist(an), ncol = length(studies), nrow = 3, byrow = F)
    dn <- lapply(1:length(studies), function(i) sub_sample_resd(network, dir, labels[[i]], gene_lists[[i]], nsub, nr,
        file))
    dn <- matrix(unlist(dn), ncol = length(studies), nrow = 3, byrow = F)

    colnames(an) <- studies
    rownames(an) <- c("res", "rand", "p")
    colnames(dn) <- studies
    rownames(dn) <- c("res", "rand", "p")

    output <- list(an, dn)
    return(output)
}
