#' Get sub-network
#'
#' @param deg A matrix (differential expression analysis)
#' @param network A matrix  (co-expression or other network)
#' @param method top (fc only) or reg (fc and q)
#' @param n_top number of genes to consider
#' @param q q-value
#' @param fct abs fold change threshold
#' @return \code{combined_output}
#' @examples
#'
#' network <- diag(1000)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:1000 )
#' colnames(network) <-paste0('gene', 1:1000 )
#' i <- 1
#' studies <- paste0('study', 1:10 )
#' nettype <- 'random'
#' dir <- 'out_test'
#' label <- studies[i]
#' filename <- 'output'
#'
#'
#' counts <- matrix(  2 ^ rpois(1e5, 5) , ncol=10, nrow=1000)
#' rownames(counts) <- paste0('gene', 1:1000 )
#' colnames(counts) <- paste0('sample', 1:10 )
#'
#' groups <- c( rep(1,5), rep(2,5) )
#' method <- 'wilcox'
#'
#' deg <- calc_DE(counts, groups, method)
#'
#' subset_network(deg$degs, network)
#'
#' @import stats utils
#' @export
#'

subset_network <- function(deg, network, method = "top",
                           n_top = 100, q = 0.05, fct = 2,
                           min_cpm=1) {

    deg <- deg[is.finite(deg$log2_fc) ,]

    m <- match(rownames(deg), rownames(network))
    f_m <- !is.na(m)
    f_am <- m[f_m]

    temp_net <- network[f_am, f_am]

    deg_sub <- deg[f_m,]

    m_cpm <- log2(deg_sub[, 1])
    fc    <- deg_sub[, 2]
    padj  <- p.adjust(deg_sub[,3], method="fdr" )

    sub_net <- list()
    deg_sig <- list()
    fc_sig <- list()
    node_degrees <- list()

    if (method == "top") {
        fc_temp = fc
        fc_temp[ m_cpm <= log2(min_cpm)  ] = 0
        deg_sig[["up"]] <- tail(order(fc_temp), n = n_top)
        deg_sig[["down"]] <- head(order(fc_temp), n = n_top)
    }
    if (method == "reg") {
        deg_sig[["down"]] <- which(fc <= -fct & padj <= q & m_cpm > log2(min_cpm) )
        deg_sig[["up"]] <- which(fc >= fct & padj <= q & m_cpm > log2(min_cpm) )
    }


    fc_sig[["down"]] <- cbind(m_cpm, fc, padj)[deg_sig[["down"]], ]
    fc_sig[["up"]]   <- cbind(m_cpm, fc, padj)[deg_sig[["up"]], ]


    sub_net[["down"]] <- temp_net[deg_sig[["down"]], deg_sig[["down"]]]
    sub_net[["up"]]   <- temp_net[deg_sig[["up"]], deg_sig[["up"]]]

    node_degrees[["down"]] <- cbind( rowSums( sub_net[["down"]]), colSums(temp_net)[deg_sig[["down"]]] )
    node_degrees[["up"]]   <- cbind( rowSums( sub_net[["up"]]), colSums(temp_net)[deg_sig[["up"]]] )
    node_degrees[["n_genes_total"]] <- dim(temp_net)[1]
    node_degrees[["n_genes_down"]]  <- length(deg_sig[["down"]])
    node_degrees[["n_genes_up"]]    <- length(deg_sig[["up"]] )


    combined_output <- list(deg_sig, fc_sig, sub_net, node_degrees)
    names(combined_output) <- c("deg_sig", "fc_sig", "sub_net", "node_degrees")
    return(combined_output)



    return(combined_output)
}
