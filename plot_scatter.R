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
