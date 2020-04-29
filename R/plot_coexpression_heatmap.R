#' Plot coexpression heatmap
#'
#' @param coexp coexpression matrix
#' @param cluster_output list
#' @param col_map color palette (array)
#' @param filt boolean to show filtered genes based on filtMin
#' @param filt_min minimum cluster size to mark for filtering
#'
#'
#' @import stats utils graphics viridis
#' @importFrom gplots heatmap.2
#' @export
#'

plot_coexpression_heatmap <- function(coexp, cluster_output,
                                      col_map = viridis(100),
                                      filt = FALSE, filt_min = 6) {


    m <- match( rownames(coexp), cluster_output$clusters[,1] )
    temp_col <- as.character(cluster_output$clusters[m,4])

    clust_size <- plyr::count( cluster_output$clusters$labels )
    clust_keep <-  clust_size[clust_size[,2] < filt_min ,1]
    genes_keep <- !is.na(match( cluster_output$clusters$labels[m], clust_keep))

    temp_filt_col <- temp_col
    if(filt==TRUE){ temp_filt_col[!genes_keep] <- "white" }

    heatmap.2(coexp, density.info = "none", trace = "none",
              col = col_map,
              Rowv = cluster_output$dendrogram,
              Colv = cluster_output$dendrogram,
              RowSideColors =  temp_col,
              ColSideColors = temp_filt_col,
              cexRow = 0.5, cexCol = 0.5, main="" )

}
