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
