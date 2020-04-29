#' Connectivity tests
#'
#' @param network A matrix (network)
#' @param dir output directory
#' @param label name
#' @param gene_list gene set of interest
#' @param flag_sub boolean, downsample
#' @return \code{output}
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)
#' studies <- paste0('study', 1:10 )
#' studies_genes   <-  paste0('gene', 1:100 )
#' colnames(genesets) <- studies
#' rownames(genesets) <- studies_genes
#' i <- 3
#' gene_list <- as.matrix(cbind(studies_genes,
#' as.numeric(genesets[,i]))[genesets[, i] > 0, ])[,1]
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:1000 )
#' colnames(network) <-paste0('gene', 1:1000 )
#' diag(network)
#' nettype <- 'random'
#' dir <- 'out_test'
#' label <- studies[i]
#' # residual_connectivity_score(network, dir, label, gene_list)
#'
#' @import stats utils graphics
#' @export

residual_connectivity_score <- function(network, dir = ".", label, gene_list, flag_sub = FALSE) {

    if (length(gene_list) < 0) {
        return(list(0, 0, 1))
    }
    n_genes <- dim(network)[1]

    m <- match(gene_list, rownames(network))
    f_r <- !is.na(m)
    f_ar <- m[f_r]


    if (sum(f_r) <= 2) {
        return(list(0, 0, 1))
    }

    network_sub <- network[f_ar, f_ar]

    node_degree <- rowSums(network, na.rm = T)
    node_degree_sub <- rowSums(network_sub, na.rm = T)

    n_genes <- dim(network)[1]
    n_genes_sub <- dim(network_sub)[1]

    if (n_genes_sub <= 2) {
        return(list(0, 0, 1))
    }

    res_sub <- residuals(node_degree[f_ar], node_degree_sub,
                         -length(node_degree_sub), n_genes, 0)

    plot(node_degree[f_ar], node_degree_sub, pch = 19,
         xlim = c(0, n_genes), ylim = c(0, n_genes_sub))
    abline(0, n_genes_sub / n_genes)
    segments(node_degree[f_ar], node_degree_sub - res_sub,
             node_degree[f_ar], node_degree_sub, col = 4)

    file <- paste0(dir, "/", paste("random", label, sep = "."),
                   "residuals.random")
    if (file.exists(file)) {
        res_rand <- unlist(read.table(file))
    } else {
        res_rand <- unlist(explore_sub_network_list_random(network,
                                dir, paste("random", label, sep = "."),
                                n_genes_sub))
    }
    res_rand <- as.matrix(res_rand)
    res_rand <- sort(res_rand)

    x <- cbind(node_degree[f_ar], node_degree_sub)
    h <- x %*% solve(t(x) %*% x) %*% t(x)
    x <- res_sub / (sd(res_sub) * sqrt(1 - diag(h)))

    a <- round(mean(res_rand, na.rm = T), 3)
    d <- round(mean(x, na.rm = T), 3)
    test <- wilcox.test(res_sub, res_rand, alt = "g")

    pvals <- sapply(1:length(x), function(i)
        sum(res_rand > x[i])) / length(res_rand)
    pvals_adj <- p.adjust(pvals, method = "BH")
    print(c(a, d, test$p.value))

    output <- cbind(as.character(rownames(network_sub)), node_degree[f_ar],
                         node_degree_sub, res_sub, x, pvals, pvals_adj)

    colnames(output) <- c("Gene", "Node degree full", "Node degree sub",
                          "Residuals", "X", "P-vals", "P-vals adj")
    write.table(output, file = paste(dir, "/", label, ".residuals", sep = ""))
    write.table(c(n_genes, n_genes_sub),
                file = paste(dir, "/", label, ".Ns", sep = ""))

    output <- list(a, d, test$p.value)
    return(output)

}
