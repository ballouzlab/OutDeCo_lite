#' Plot network
#'
#' @param sub_net network of genes
#' @param clust_net output from cluster_coexp
#' @param threshold weighted edge threshold to plot
#' @param filt_min minimum cluster size to display in plot
#' @examples
#'
#' network <- diag(100)
#' upper <- row(network) < col(network)
#' network[upper] <- runif(sum(upper), 0,1 )
#' network <- network + t(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:100 )
#' colnames(network) <-paste0('gene', 1:100 )
#' clust_net <- cluster_coexp(network)
#' plot_network(network, clust_net)
#'
#' @import viridis venn 
#' @importFrom gplots heatmap.2
#' @importFrom igraph graph_from_data_frame V E delete_edges
#' @export

plot_network <- function(sub_net, clust_net, threshold = 0.5, filt_min = 6) {

  diag(sub_net) <-  0
  upper <- row(sub_net) < col(sub_net)
  pairs <- which(upper, arr.ind = T )
  gene_names <- rownames(sub_net)
  weights <- sub_net[upper]
  pairs <- data.frame( p1 = gene_names[pairs[,1]], p2 = gene_names[pairs[,2]] , weights = weights )


  # inet <- igraph::graph_from_adjacency_matrix(sub_net, weighted = T, mode = "undirected")
  inet <- igraph::graph_from_data_frame(pairs, directed=F)
  igraph::E(inet)$weight <- 1 - weights
  igraph::E(inet)$width <- (weights^2 * 10)
  igraph::E(inet)$edge.color <- "gray80"
  igraph::E(inet)$color <- viridis(100)[ round(weights * 100) ]

    o <- match(igraph::V(inet)$name, clust_net$clusters$genes)

  igraph::V(inet)$color <- as.character(clust_net$clusters$colors)[o]

  igraph::V(inet)$label <- ""
  igraph::V(inet)$size <- 4

  clust_size <- plyr::count(clust_net$clusters$labels )
  clust_keep <-  clust_size[clust_size[,2] < filt_min ,1]
  genes_keep <- !is.na(match( clust_net$clusters$labels, clust_keep))

  o <- match( igraph::V(inet)$name, clust_net$clusters$genes[genes_keep])
  f_n <- !is.na(o)
  igraph::V(inet)$size[f_n] <- 10


  inet_sub <-  igraph::delete_edges(inet, igraph::E(inet)[weights < threshold])

  plot(inet_sub )
  #, layout = layout_with_fr )
  return(inet_sub)
}




