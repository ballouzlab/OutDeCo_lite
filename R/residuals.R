#' residuals from identity line
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param a maximum x value
#' @param b maximum y value
#' @param c y intercept
#' @return \code{res}
#'
#' @examples
#' x <- (rnorm(1000))
#' y <- (jitter(x)) + (rnorm(1000))
#' rs <- residuals(x,y, -1,1,0)
#'
#' plot(x,y, pch=19, cex=1.5)
#' abline(0,1, col='grey', lty=2, lwd=2)
#' segments(x, y - rs, x ,y , col=4)
#'
#' @export
#'


residuals <- function(x, y, a, b, c) {
    res <- (a * x + b * y + c)
    return(res)
}
