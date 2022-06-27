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
                                      filt = FALSE, filt_min = 6, 
                                      flag_plot_bin= FALSE) {


    m <- match( rownames(coexp), cluster_output$clusters[,1] )
    temp_col <- as.character(cluster_output$clusters[m,4])

    clust_size <- plyr::count( cluster_output$clusters$labels )
    clust_keep <-  clust_size[clust_size[,2] < filt_min ,1]
    genes_keep <- !is.na(match( cluster_output$clusters$labels[m], clust_keep))

    temp_filt_col <- temp_col
    if(filt==TRUE){ temp_filt_col[!genes_keep] <- "white" }
    
    if(flag_plot_bin == FALSE) { 
        temp = coexp
    } else {
       temp = cluster_output$distance_matrix
    } 
     
     heatmap.2(temp, density.info = "none", trace = "none",
              col = col_map,
              Rowv = cluster_output$dendrogram,
              Colv = cluster_output$dendrogram,
              RowSideColors =  temp_col,
              ColSideColors = temp_filt_col,
              cexRow = 0.5, cexCol = 0.5, main="" )  
  
}




#' Plot gene set enrichment results
#'
#' @param data gene set enrichment results output
#' @param gene_list list of genes used as input into the gene set enrichment
#' @param gene_sets gene sets used as input into the gene set enrichment analysis
#'
#' @import stats utils graphics viridis RColorBrewer
#' @importFrom gplots heatmap.2
#' @export
#'

plot_gene_set_enrichment <- function(data, gene_list, gene_sets) {

   colssig = gplots::colorpanel(100, "lightgrey", "red", "darkmagenta")

   paths <- data$term
   paths_padj <- data$padj
   paths_pval <- data$pval

   # Match gene list with genes in gene sets
   m <-  match( gene_list, rownames(gene_sets) )
   sub_gene_sets <- gene_sets[m,]
   sub_gene_sets[is.na(sub_gene_sets)] = 0
   rownames(sub_gene_sets) <- gene_list

   # Match enrichment results with gene sets input
   m <- match( paths, colnames(gene_sets) )
   sub_gene_sets <- sub_gene_sets[,m]

   # Set up pval and padj matrices to plot
   sub_gene_sets_pval <- (sub_gene_sets * 0 ) + 1
   sub_gene_sets_padj <- (sub_gene_sets * 0 ) + 1

   for(i in 1:length(paths_pval)){
       sub_gene_sets_pval[ sub_gene_sets[,i]==1,i] = paths_pval[i]
       sub_gene_sets_padj[ sub_gene_sets[,i]==1,i] = paths_padj[i]
   }

   log10sub_gene_sets_pval = -log10(sub_gene_sets_pval)
   log10sub_gene_sets_padj = -log10(sub_gene_sets_padj)
   log10sub_gene_sets_pval[!is.finite(log10sub_gene_sets_pval)] = 0
   log10sub_gene_sets_padj[!is.finite(log10sub_gene_sets_padj)] = 0

  filt = colSums(log10sub_gene_sets_pval) > 0
 #  gplots::heatmap.2(log10sub_gene_sets_pval[,]  , Rowv=F, Colv=F,

  gplots::heatmap.2(log10sub_gene_sets_pval[,filt]  , Rowv=F, Colv=F,
              col=colssig, cexRow = 0.7, cexCol = 0.7,
              notecol="black",
              notecex=1,
              keysize=1,
              key.xlab="-log10 adjusted P-value",
              key.title="Enrichment", trace="none", density="none" )

}


plot_gene_set_enrichment_ranked <- function(data, gene_rankings, gene_list, gene_sets, psig = 0.05) {

  par(mfrow=c(2,2))
  n <- dim(gene_sets)[1]
  nn <- dim(gene_sets)[2]

  rocs <- lapply(1:nn, function(i) get_roc( gene_rankings, gene_sets[, i]))
  pvals <- unlist(lapply(1:nn, function(i) wilcox.test( gene_rankings[gene_sets[, i]==1],  gene_rankings[gene_sets[, i]==0])$p.val ))
  # pvals <- unlist(lapply(1:nn, function(i) wilcox.test( gene_rankings[gene_sets[, i]==1],  gene_rankings[gene_sets[, i]==0], alt="g")$p.val ))


  padj <- p.adjust(pvals)
  sig = (padj<psig) * 1  + 1
  names(sig) = colnames(gene_sets)


  # Panel 1
  plot(-10,-10, xlim=c(0,1), ylim=c(0,1),xlab = "FPR", ylab = "TPR", bty = "n"   )
  ll <- lapply(1:nn, function(i) lines( rocs[[i]][,1],rocs[[i]][,2] , col=EGAD::make_transparent(sig[i],100), lwd=sig[i] ))
  ll <- lapply(1:which(sig==2), function(i) lines( rocs[[i]][,1],rocs[[i]][,2] , col=EGAD::make_transparent(sig[i],100), lwd=sig[i] ))

  # Panel 2
  nsets  <- colSums(gene_sets)

  hist(nsets, breaks=30, col=viridis(10)[5], border=NA, xlab="Gene set size", main="")

  # Panel 3
  hist(data, breaks=30, col=viridis(10)[5], border=NA, xlab="Gene set AUCs", main="")

  # Panel 4
  plot( data, -log10(pvals), pch=19, bty="n", xlab="Gene set AUCs", ylab="-log10 P-value", cex=sig, col=sig)


  # print( cbind( sig, pvals, padj, nsets))

}





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



#' Plot scatterplot and histograms
#'
#' @param x values
#' @param y values
#' @param clusters clusters if applicable
#' @param xlab x axis labe
#' @param ylab y axis label
#' @param main main title
#' @param xybreaks number of breaks for histogram
#' @param flag type of plot to draw (hist or density)
#' @param ... plotting parameters
#' @examples
#'#'
#' x <- (rnorm(1000))
#' y <- (jitter(x)) + (rnorm(1000))
#' plot_scatter_hist(x,y)
#'
#' @import viridis stats utils graphics
#' @importFrom scales alpha
#' @export

plot_scatter <- function(x, y, clusters = FALSE, xlab = "", ylab = "", main = "",
                                 xybreaks = 100,  flag ="hist", ...) {

    par_old <- par()

    if( flag == "density") {
        plot_scatter_density(x, y, clusters, xlab = xlab, ylab = ylab, main = main ,
                             xybreaks = xybreaks,  ...)
    }
    else {
       plot_scatter_hist(x, y, clusters, xlab = xlab, ylab = ylab, main = main,
                         xybreaks = xybreaks,   ... )

    }
  par(par_old)
}





#' Plot scatterplot and histograms
#'
#' @param x values
#' @param y values
#' @param clusters output from cluster_coexp
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main main title
#' @param xybreaks number of breaks for histogram
#' @param ... plotting parameters
#' @examples
#'#'
#' x <- (rnorm(1000))
#' y <- (jitter(x)) + (rnorm(1000))
#' plot_scatter_hist(x,y)
#'
#' @import viridis stats utils graphics
#' @importFrom scales alpha
#' @export

plot_scatter_density <- function(x, y, clusters = FALSE, clust_flag=FALSE, xlab = "", ylab = "", main = "",
                                 xybreaks = 100,  ...) {

    par(oma = c(3, 3, 0, 0))
    xyrange <- c(floor(min(c(x, y))), ceiling(max(c(x, y))))
    zones <- matrix(c(2, 0, 1, 3), ncol = 2, byrow = TRUE)
    layout(zones, widths = c(4/5, 1/5), heights = c(1/5, 4/5))
    x_max <- max(x)
    y_max <- max(y)
    x_min <- min(x)
    y_min <- min(y)


    if( length(clusters) > 1 ){
        n_clusters <- max(clusters$labels)
        col_points <- scales::alpha( as.character(clusters$colors) ,0.8)
        color_legend <-  unique(cbind(clusters$labels,as.character(clusters$colors)))
        score_x_list <- lapply(1:n_clusters, function(j) x[clusters$labels==j])
        score_y_list <- lapply(1:n_clusters, function(j) y[clusters$labels==j])
    }

    else {
      n_clusters <- 1
      col_points <- scales::alpha("black",0.5)
      color_legend <-  cbind (1, "grey")
      score_x_list <- lapply(1:n_clusters, function(j) x)
      score_y_list <- lapply(1:n_clusters, function(j) y)


    }


    x_dens <- lapply(n_clusters:1, function(i) density(score_x_list[[i]], bw=0.01,  from=x_min, to=x_max)  )
    y_dens <- lapply(n_clusters:1, function(i) density(score_y_list[[i]], bw=0.01,  from=y_min, to=y_max)  )
    x_dens_sums <- as.matrix( sapply(n_clusters:1, function(i) x_dens[[i]]$y ) )
    y_dens_sums <- as.matrix( sapply(n_clusters:1, function(i) y_dens[[i]]$y ))
  #  dens_max = max( c(rowSums(x_dens_sums), rowSums(y_dens_sums)))


    # Plot main scatter
    par(mar = c(3, 3, 1, 1))
    plot(x, y, pch = 21, main = main, xlim=c(x_min,x_max), ylim=c(y_min,y_max),
         axes=F, col= col_points, lwd=2, bg = "lightgray", ...)
     abline(0, 1, col = "grey", lwd = 3)
    axis(1)
    axis(2)
    abline( v = mean(x, na.rm=T), col = "grey", lty=2)
    abline( h = mean(y, na.rm=T), col = "grey", lty=2)


    # Plot x-axis density
    par(mar = c(0, 3, 1, 1))

    dens_max = max(  rowSums(x_dens_sums) )
    plot(x_dens[[1]]$x, x_dens[[1]]$y,
         xlab="", ylab="",
         xlim=c(x_min,x_max), ylim=c(0,dens_max),
         col=0, axes=F)
    axis(2)
    if (n_clusters > 1) {
      for (i in n_clusters:2) {
        polygon(
          c(x_min, x_dens[[1]]$x, x_max),
          c(0, rowSums(x_dens_sums[, 1:i]), 0) ,
          col = (color_legend[, 2])[i],
          border = NA
        )
      }

    }
    polygon( c(x_min,x_dens[[1]]$x,x_max),
             c(0,x_dens_sums[,1],0) ,
             col=(color_legend[,2])[1], border=NA)


    # Plot y-axis density
    par(mar = c(3, 0, 1, 1))

    dens_max = max(  rowSums(y_dens_sums) )
     plot(  y_dens[[1]]$y, y_dens[[1]]$y,
           xlab="", ylab="",
           ylim=c(y_min,y_max), xlim=c(0,dens_max),
           col=0, axes=F)
    axis(1)

    if (n_clusters > 1) {
      for (i in n_clusters:2) {
        polygon(
          c(0, rowSums(y_dens_sums[, 1:i]), 0) ,
          c(y_min, y_dens[[1]]$x, y_max),
          col = (color_legend[, 2])[i],
          border = NA
        )
      }
    }
    polygon( c(0,y_dens_sums[,1],0) ,
             c(y_min,y_dens[[1]]$x,y_max),
             col=(color_legend[,2])[1], border=NA)

    # Add x and y axes
    par(oma = c(3, 3, 0, 0))
    mtext(xlab, side = 1, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(x) - min(x))/(max(x) - min(x)))
    mtext(ylab, side = 2, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(y) - min(y))/(max(y) - min(y)))
}




#' Plot scatterplot and histograms
#'
#' @param x values
#' @param y values
#' @param xlab x axis labe
#' @param ylab y axis label
#' @param main main title
#' @param xybreaks number of breaks for histogram
#' @param clusters output from cluster_coexp
#' @param ... plotting parameters
#' @examples
#'#'
#' x <- (rnorm(1000))
#' y <- (jitter(x)) + (rnorm(1000))
#' plot_scatter_hist(x,y)
#'
#' @import viridis stats utils graphics
#' @importFrom scales alpha
#' @export

plot_scatter_hist <- function(x, y, xlab = "", ylab = "", main = "", xybreaks = 100,
                              clusters = FALSE, ...) {

    par(oma = c(3, 3, 0, 0))
    xyrange <- c(floor(min(c(x, y))), ceiling(max(c(x, y))))
    zones <- matrix(c(2, 0, 1, 3), ncol = 2, byrow = TRUE)
    layout(zones, widths = c(4/5, 1/5), heights = c(1/5, 4/5))
    xhist <- hist(x, plot = FALSE, breaks = xybreaks)
    yhist <- hist(y, plot = FALSE, breaks = xybreaks)

    top <- max(c(xhist$counts, yhist$counts))

    if(   length(clusters)>1 ) {
        col_points <- scales::alpha( as.character(clusters$colors) ,0.8)
        n_clusters <- max(clusters$labels)
        color_legend <- unique(cbind(clusters$labels,as.character(clusters$colors)))
        xhist_mat <- sapply(1:n_clusters, function(i) hist(x[clusters$labels == i], plot = FALSE, breaks = xhist$breaks)$counts )
        yhist_mat <- sapply(1:n_clusters, function(i) hist(y[clusters$labels == i], plot = FALSE, breaks = yhist$breaks)$counts )

    }
    else {
        col_points <- scales::alpha("black",0.5)
        }


    # Plot main scatter
    par(mar = c(3, 3, 1, 1))
    plot(x, y, pch = 21, main = main, axes=F, col= col_points, lwd=2, bg = "lightgray", ...)
    abline(0, 1, col = "grey", lwd = 3)
    axis(1)
    axis(2)
    abline( v = mean(x, na.rm=T), col = "grey", lty=2)
    abline( h = mean(y, na.rm=T), col = "grey", lty=2)

    # Plot x-axis histogram
    par(mar = c(0, 3, 1, 1))
    if(  length(clusters)>1 ) {
        bp <- barplot( t(xhist_mat), axes= FALSE, ylim = c(0, top), space = 0,
                       col= color_legend[,2], border=NA )
        }
    else { bp <- barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,
                         col = "grey", border=NA)
    }

    axis(2)
    x1 <- xhist$mids[1]
    y1 <- bp[1]
    m  <- diff(bp)[1] / diff(xhist$mids)[1]
    b  <- y1 - m * x1


    abline( v = m*mean(x, na.rm=T)+b, col = "grey", lty=2)

    # Plot y-axis histogram
    par(mar = c(3, 0, 1, 1))

    if(  length(clusters)>1 ) {
        bp <- barplot( t(yhist_mat), axes= FALSE, xlim = c(0, top), space = 0,
                       col= color_legend[,2],horiz = TRUE, border=NA )
    }
    else { bp <- barplot( yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,
                         col = "grey", horiz = TRUE,border=NA)
    }


    axis(1)
    x1 <- yhist$mids[1]
    y1 <- bp[1]
    m  <- diff(bp)[1] / diff(yhist$mids)[1]
    b  <- y1 - m * x1

    abline( h = m * mean(y, na.rm=T) + b , col = "grey", lty=2)

    # Add x and y labels
    par(oma = c(3, 3, 0, 0))
    mtext(xlab, side = 1, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(x) - min(x))/(max(x) - min(x)))
    mtext(ylab, side = 2, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(y) - min(y))/(max(y) - min(y)))

}


