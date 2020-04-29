#' Plot 2D histogram
#'
#' @param mat yA matrix of frequencies
#' @param pt.x threshold x
#' @param pt.y threshold
#' @param ... parameters to pass onto the image() function
#' @examples
#'
#' data_pre <- as.matrix(rowSums(matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)))
#' data_post <- data_pre - (data_pre * rbinom(100, 1, 0.2 ))
#' recurs = cbind(data_post, data_pre)
#'
#' mat10 <- get_recur_mat(recurs)
#' plot_2D_hist(mat10, 1, 4)
#'
#' @import stats utils graphics
#' @export
#'

plot_2D_hist <- function(mat, pt.x, pt.y, ...) {
    image(mat, axes = F, ...)
    n <- dim(mat)[2] - 1
    ni <- diff((0:n/n))[1]
    abline(h = (0:n/n)[pt.y] + ni/2, lwd = 3, col = "grey")
    axis(2, at = 0:n/n, labels = 0:n)
    n <- dim(mat)[1] - 1
    ni <- diff((0:n/n))[1]
    axis(1, at = 0:n/n, labels = (0:n))
    abline(v = (0:n/n)[pt.x] + ni/2, lwd = 3, col = "grey")
}
