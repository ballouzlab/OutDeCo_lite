#' Plot recurrence histogram
#'
#' @param recur values
#' @param fdrs output from calc_fdrs_recur
#' @param n_max total number of studies
#' @param flag_plot plot type
#' @examples
#'
#' data <-  matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' n_max <- dim(data)[2]
#' fdrs <- calc_fdrs_recur( data )
#' recur <- rowSums(data, na.rm=TRUE)
#' plot_recurrence( data, fdrs, n_max)
#'
#' @import viridis venn
#' @importFrom gplots colorpanel
#' @export

plot_multi_matrix <- function(data, fdrs, n_max, flag_plot = "hist", ...) {

   recur <- rowSums(data, na.rm=TRUE)
   data <- data[recur>0,]


    col_map <- gplots::colorpanel(n_max+1, "white",  "darkmagenta")
    o_c <- order (colSums(data))
    data2 <- cbind(rowSums(data), data)
    o_r <- rev(do.call( order, lapply(1:(n_max+1), function (i) data2[,i]) ) )
    n_genes  <- dim(data)[1]


    par(oma = c(1, 7, 1, 1))
    par(mar = c(1, 7, 1, 1))
    zones <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
    layout(zones, widths = c(1,1,1,1), heights = c(1, 1, 1,1))

    image( as.matrix(data2[o_r,1]),  col = 0, axes=F )
    image( as.matrix(data2[o_r,1]),  col = col_map, axes=F )
    mtext( "Recurrence", side = 2, at = 0.25 , las=1, cex=1.5 )

    par(mar = c(0, 7, 0, 1))
    image( (data[o_r,o_c]),  col = c(0,1), axes=F )
    mtext( colnames(data)[o_c], side = 2, at =  0:(n_max-1 )  / (n_max-1) , las=1, cex=0.75 )



}




