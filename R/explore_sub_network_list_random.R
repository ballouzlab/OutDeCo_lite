#' Connectivity tests
#'
#' @param network A matrix (network)
#' @param dir output directory
#' @param label name
#' @param ns number of genes to subsample to
#' @param nr number of repeats
#' @return \code{data}
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)
#' studies <- paste0('study', 1:10 )
#' studies_genes   <-  paste0('gene', 1:100 )
#' colnames(genesets) <- studies
#' rownames(genesets) <- studies_genes
#' i <- 1
#' gene_list <- as.matrix(cbind(studies_genes, as.numeric(genesets[,i]))[genesets[, i] > 0, ])[,1]
#' ns <- length(gene_list)
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:1000 )
#' colnames(network) <-paste0('gene', 1:1000 )
#' nettype <- 'random'
#' dir <- '.'
#'
#' # explore_sub_network_list_random(network, 'output', 'test', ns)
#'
#'
#' @import stats utils
#' @export
#'
#'
#'

explore_sub_network_list_random <- function(network, dir, label, ns, nr = 1000) {
    
    node_degree = rowSums(network, na.rm = T)
    # print( mean(node_degree, na.rm=T) ) print ( sd(node_degree, na.rm=T) )
    N = dim(network)[1]
    summary = matrix(0, ncol = 2, nrow = nr)
    data = list()
    
    for (i in 1:nr) {
        
        gene_list = sample(rownames(network), ns)
        m = match(gene_list, rownames(network))
        f.r = !is.na(m)
        f.ar = m[f.r]
        
        network_sub = network[f.ar, f.ar]
        node_degree_sub = rowSums(network_sub, na.rm = T)
        
        n = dim(network_sub)[1]
        x = c(0, N)
        y = c(n/N * x[1], n)
        k = 1
        data.sub = list()
        frac = 1
        prb = sample(gene_list[f.r], length(gene_list[f.r])/frac)
        m = match(gene_list[f.r], prb)
        prb = !is.na(m)
        network_sub.prb = network_sub[prb, prb]
        node_degree_sub.prb = rowSums(network_sub.prb, na.rm = T)
        X = residuals(node_degree[f.ar][prb], node_degree_sub.prb, -length(node_degree_sub.prb), N, 0)
        H = X %*% solve(t(X) %*% X) %*% t(X)
        resd.sub.prb = X/(sd(X) * sqrt(1 - diag(H)))
        summary[i, k] = mean(resd.sub.prb, na.rm = T)
        summary[i, k + 1] = sd(resd.sub.prb, na.rm = T)
        
        data[[i]] = resd.sub.prb
        # data[[i]] = X
        
    }
    write.table(summary, file = paste(dir, "/", label, "summary.stats.random", sep = ""))
    write.table(matrix(unlist(data), nrow = nr, byrow = T), file = paste(dir, "/", label, "residuals.random", sep = ""))
    return(data)
}
