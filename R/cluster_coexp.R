#' Cluster sub-network
#'
#' @param coexp A matrix (co-expression or other network)
#' @param runid ID title
#' @param filt_min minimum size of cluster to consider
#' @param medK median expression value of co-expression network
#' @param flag_plot Boolean
#' @param col_map color palette for heatmaps
#' @return \code{clusterids}
#' @examples
#'
#' cpm <- matrix(  abs(rnorm(100)*10) , ncol=1, nrow=100)
#' cpm <- sapply(1:10, function(i) cpm[,1])
#' cpm <- jitter(cpm, amount = 4)
#' network <- matrix(rank(cor(t(cpm))), nrow=100) - 1
#' network <- network/max(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:100 )
#' colnames(network) <-paste0('gene', 1:100 )
#' nettype <- 'random'
#' dir <- 'out_test'
#' cluster_coexp(network)
#'
#' @import dynamicTreeCut viridis stats utils
#' @importFrom gplots heatmap.2
#' @export
#'
#'

cluster_coexp <- function(coexp, runid = "", filt_min = 6, medK = 0.5,
                                 col_map = viridis(100), method = "average",
                                 flag_plot = FALSE, flag_med = TRUE, flag_dist = FALSE,
                                 frac = 0.995, deep_split = 2, min_cs = 2 ) {

    # Get distance matrix from co-expression
    temp <- coexp

    # Make as a binary network based on median
    if( flag_med == TRUE) {
        temp[temp > medK] <- 1
        temp[temp <= medK] <- 0
        diag(temp) <- 0
    }

    gene_names <- rownames(temp)

    # Re-calculate distances between genes for distance matrix
    if( flag_dist == TRUE) {
        dist_temp <- dist(temp)
    } else {
        dist_temp <- as.dist(temp)
    }

    # Cluster genes using distance matrix
    clust_tree <- hclust(dist_temp, method = method)
    clust_dend <- as.dendrogram(clust_tree)

    if( flag_dist == TRUE ){
        max_h = max(clust_tree$height)
    } else {
        max_h = 1
    }

    # Extract clusters/modules
    unmerged_modules <- dynamicTreeCut::cutreeDynamic(dendro = clust_tree,
                                                      distM = as.matrix(dist_temp),
                                                      deepSplit = deep_split,
                                                      cutHeight = frac * max_h,
                                                      minClusterSize = min_cs,
                                                      pamRespectsDendro = FALSE)

    n_max <- max(unmerged_modules)
    n_l <- sum(unmerged_modules==0)
    merged_modules <- unmerged_modules
    if( n_l > 0) { merged_modules[unmerged_modules==0] <- (1:n_l ) + n_max }
    merged_modules <- merged_modules[clust_tree$order]

    # Re-label clusters ids
    i.prev <- ""
    ki <- 1
    ji <- 0
    remerged_modules <- as.numeric(merged_modules) * 0
    for (ii in as.numeric(merged_modules) ) {
        if (ii == i.prev) {
            remerged_modules[ki] <- ji
        } else {
            i.prev <- ii
            ji <- ji + 1
            remerged_modules[ki] <- ji

        }
        ki <- ki + 1
    }

    # Total number of clusters
    nsclust <- as.numeric(remerged_modules) + 1


    # Generate colors for modules
    merged_colors <- viridis::magma(max(nsclust))[nsclust]
    m <- match( gene_names, gene_names[clust_tree$order] )

    # Plot
    if (flag_plot == TRUE) {

        heatmap.2( coexp, density.info = "none", trace = "none",
                  col = col_map,
                  Rowv = clust_dend, Colv = clust_dend,
                  RowSideColors = merged_colors[m],
                  ColSideColors = merged_colors[m],
                  cexRow = 0.5, cexCol = 0.5, main = runid)


    }

    # Tidy up output
    clusterids <- data.frame( genes = gene_names[clust_tree$order], labels = remerged_modules,
                             labels_unmerged = unmerged_modules[clust_tree$order],
                             colors = merged_colors)

    combined_output <- list(as.matrix(dist_temp), clust_tree, clust_dend, m, clusterids)
    names(combined_output) <- c( "distance_matrix", "tree", "dendrogram", "order", "clusters")
    return(combined_output)
}
