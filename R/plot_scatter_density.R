#' Plot scatterplot and histograms
#'
#' @param x values
#' @param y values
#' @param xlab x axis labe
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

plot_scatter_density <- function(x, y, clusters = FALSE, xlab = "", ylab = "", main = "",
                                 xybreaks = 100,  ...) {

    par(oma = c(3, 3, 0, 0))
    xyrange <- c(floor(min(c(x, y))), ceiling(max(c(x, y))))
    zones <- matrix(c(2, 0, 1, 3), ncol = 2, byrow = TRUE)
    layout(zones, widths = c(4/5, 1/5), heights = c(1/5, 4/5))
    x_max <- max(x)
    y_max <- max(y)

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


    x_dens <- lapply(n_clusters:1, function(i) density(score_x_list[[i]],   from=0, to=x_max)  )
    y_dens <- lapply(n_clusters:1, function(i) density(score_y_list[[i]],   from=0, to=y_max)  )
    x_dens_sums <- as.matrix( sapply(n_clusters:1, function(i) x_dens[[i]]$y ) )
    y_dens_sums <- as.matrix( sapply(n_clusters:1, function(i) y_dens[[i]]$y ))
  #  dens_max = max( c(rowSums(x_dens_sums), rowSums(y_dens_sums)))


    # Plot main scatter
    par(mar = c(3, 3, 1, 1))
    plot(x, y, pch = 21, main = main, xlim=c(0,x_max), ylim=c(0,y_max),
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
         xlim=c(0,x_max), ylim=c(0,dens_max),
         col=0, axes=F)
    axis(2)
    if (n_clusters > 1) {
      for (i in n_clusters:2) {
        polygon(
          c(0, x_dens[[1]]$x, x_max),
          c(0, rowSums(x_dens_sums[, 1:i]), 0) ,
          col = (color_legend[, 2])[i],
          border = NA
        )
      }

    }
    polygon( c(0,x_dens[[1]]$x,x_max),
             c(0,x_dens_sums[,1],0) ,
             col=(color_legend[,2])[1], border=NA)


    # Plot y-axis density
    par(mar = c(3, 0, 1, 1))

    dens_max = max(  rowSums(y_dens_sums) )
     plot(  y_dens[[1]]$y, y_dens[[1]]$y,
           xlab="", ylab="",
           ylim=c(0,y_max), xlim=c(0,dens_max),
           col=0, axes=F)
    axis(1)

    if (n_clusters > 1) {
      for (i in n_clusters:2) {
        polygon(
          c(0, rowSums(y_dens_sums[, 1:i]), 0) ,
          c(0, y_dens[[1]]$x, y_max),
          col = (color_legend[, 2])[i],
          border = NA
        )
      }
    }
    polygon( c(0,y_dens_sums[,1],0) ,
             c(0,y_dens[[1]]$x,y_max),
             col=(color_legend[,2])[1], border=NA)

    # Add x and y axes
    par(oma = c(3, 3, 0, 0))
    mtext(xlab, side = 1, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(x) - min(x))/(max(x) - min(x)))
    mtext(ylab, side = 2, line = 1, outer = TRUE, adj = 0, at = 0.5 * (mean(y) - min(y))/(max(y) - min(y)))
}
