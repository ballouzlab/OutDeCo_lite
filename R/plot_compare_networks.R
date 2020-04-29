#' A plot of two networks (and their joining river)
#'
#' @param net1 first network
#' @param net2 second network
#' @param clust_net1 clusters of first
#' @param clust_net2 clusters of second
#'
#'
#' @import viridis riverplot
#' @importFrom gplots colorpanel
#' @export

plot_compare_networks <- function(net1, net2,
                                  clust_net1, clust_net2) {
    par_old <- par()
    # Subset networks  to common genes, in case not done prior
    m =   match( rownames(net1), rownames(net2)  )
    f.n1 = !is.na(m)
    f.n2 = m[f.n1]

    if( sum(f.n1) != length(f.n1) ) {
        print("Networks aren't compatible, subset and recluster!")
        return(0)
    }
    # Reorder for heatmaps
    m1 <- match(clust_net1$clusters$genes, rownames(net1)  )
    m2 <- match(clust_net2$clusters$genes, rownames(net2)  )


    xn <- 0
    # Set up riverplot
    ## Common order between clusters for riverplot
    m = match(clust_net1$clusters[,1], clust_net2$clusters[,1] )
    f.g = !is.na(m)
    f.n = m[f.g]

    a <- clust_net1$clusters[f.g,2]
    b <- clust_net2$clusters[f.n,2] + max(a)

    a_c <- (unique(clust_net1$clusters[f.g,c(2,4)] ) )
    b_c <- (unique(clust_net2$clusters[f.n,c(2,4)] ) )
    b_c[,1] = b_c[,1] + max(a)

    b2 <- clust_net2$clusters[,2]

    temp <- rbind(a_c, b_c)


    nodes <- data.frame( ID=temp[,1],
                         x = c(rep(1, max(a)), rep(2, max(b)-max(a) )),
                         col= as.character(temp[,2]), stringsAsFactors = FALSE)
    o     <- order(nodes[,1])
    nodes <- nodes[o,]
    temp  <- data.frame( N1= a, N2= b, Value = 1 )
    freq  <- plyr::count(temp)
    edges <- freq[,1:3]
    edges[,3] <- freq[,4]
    r <- riverplot::makeRiver( nodes, edges,  node_labels="" )


    # color the same between the two networks
     nmax = round(max( net1, net2) , 2)
     nnn =  (nmax*100)
     nnn0  =  floor( min(net1,net2)*100)
     b3 = c(nnn0:nnn)/100
     hb = hist(b3, breaks=5, plot=F)

    #nmax <- 1
    #nnn  <- nmax*100
    #b3   <- c(0:nnn)/100
    #hb   <- hist(b3, breaks=5, plot=F)
    cols_net <- viridis::viridis(length(b3)-1 )

    par(oma = c(2, 2, 2, 2))
    par(mar = c(1, 1, 1, 1))
    zones <- matrix(c(1:12), ncol = 3, nrow=4, byrow = TRUE)

    layout(zones, heights = c(0.5,0.5,5, 0.5))

    par(mar = c(xn, xn, xn, xn))
    plot(clust_net1$dendrogram)
    plot(0, col=0, axes=F)
    plot(clust_net2$dendrogram)

    par(mar = c(1, 1, 1, 1))
    image(cbind(a,a), col=viridis::magma(max(a)+1)[-1], axes=F)
    plot(0, col=0, axes=F)
    image(cbind(b2,b2), col=viridis::magma(max(b2)+1)[-1], axes=F)

    par(mar = c(0, 1, 0, 1))
    image( net1[m1,m1], col = cols_net,  axes=F, breaks=b3)
    par(mar = c(0, 0, 0, 0))
    rr <- capture.output(plot( r, node_margin = 0, plot_area=1))
    par(mar = c(0, 1, 0, 1))
    image( net2[m2,m2],  col = cols_net, axes=F, breaks=b3 )

    par(mar = c(1, 5, 1, 5))
    image( cbind(b3,b3), col=cols_net, breaks=b3, axes=F  )
    axis(1, at= 0:(length(hb$breaks)-1)/(length(hb$breaks)-1) , lab=c(hb$breaks)  )

    par(par_old)
}




