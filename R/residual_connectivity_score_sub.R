#' Connectivity tests
#'
#' @param network A matrix (network)
#' @param dir output directory
#' @param label name
#' @param gene_list gene set of interest
#' @param n_sub number of genes to subsample
#' @param n_r number of repeats
#' @param fileout filename
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)
#' studies <- paste0('study', 1:10 )
#' studies_genes   <-  paste0('gene', 1:100 )
#' colnames(genesets) <- studies
#' rownames(genesets) <- studies_genes
#' i <- 3
#' gene_list <- as.matrix(cbind(studies_genes,
#'                as.numeric(genesets[,i]))[genesets[, i] > 0, ])[,1]
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <- paste0('gene', 1:1000 )
#' colnames(network) <- paste0('gene', 1:1000 )
#' nettype <- 'random'
#' dir <- 'out_test'
#' label <- studies[i]
#' fileout <- 'output'
#'
#' #residual_connectivity_score_sub(network, dir, label,
#' #gene_list, 20, 100, fileout)
#'
#' @import stats utils graphics
#' @export
#'

residual_connectivity_score_sub <- function(network, dir, label,
                                            gene_list, n_sub, n_r, fileout) {
    if (length(gene_list) < 0) {
        return(list(0, 0, 1))
    }

    m <- match(gene_list, rownames(network))
    f_r <- !is.na(m)
    f_ar <- m[f_r]

    if (sum(f_r) <= n_sub) {
        return(list(0, 0, 1))
    }

    network_sub <- network[f_ar, f_ar]
    node_degree <- rowSums(network, na.rm = TRUE)
    node_degree_sub <- rowSums(network_sub, na.rm = TRUE)
    n_genes <- dim(network)[1]
    n_genes_sub <- dim(network_sub)[1]

    if (n_genes_sub <= n_sub) {
        return(list(0, 0, 1))
    }

    if (file.exists(fileout)) {
        res_rand <- unlist(read.table(fileout))
    } else {
        res_rand <- unlist(explore_sub_network_list_random(network, dir,
                                                           paste("random.sub.",
                                                            label, sep = "."),
                                                           n_sub))
    }

    res_rand <- as.matrix(res_rand)
    res_rand <- sort(res_rand)

    mat <- matrix(0, ncol = 5, nrow = n_r)
    for (k in 1:n_r) {
        f_samp <- sample(n_genes, n_sub, replace = FALSE)
        node_degree_sub <- rowSums(network_sub[f_samp, f_samp])

        res_sub <- residuals(node_degree[f_ar][f_samp], node_degree_sub,
                             -length(node_degree_sub), n_genes, 0)

        # studentized residuals
        y <- cbind(node_degree[f_ar][f_samp], node_degree_sub)
        h <- y %*% solve(t(y) %*% y) %*% t(y)
        x <- res_sub / (sd(res_sub) * sqrt(1 - diag(h)))

        pvals <- sapply(1:length(x), function(i)
            sum(res_rand > x[i])) / length(res_rand)
        pvals_adj <- p.adjust(pvals, method = "BH")
        a <- round(mean(res_sub, na.rm = TRUE), 3)
        b <- mean(pvals, na.rm = TRUE)
        c <- mean(pvals_adj, na.rm = TRUE)
        d <- round(mean(x, na.rm = TRUE), 3)
        test <- wilcox.test(res_sub, res_rand, alt = "g")
        mat[k, ] <- c(a, b, c, d, test$p.value)


    }
    write.table(mat, file = paste(dir, label, ".resamp.residuals", sep = ""))

    a <- round(mean(res_rand, na.rm = TRUE), 3)
    d <- round(mean(mat[, 4], na.rm = TRUE), 3)

    test$p_value <- geo_mean(mat[, 5])
    output <- list(a, d, test$p_value)
    return(output)

}
