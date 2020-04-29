#' geometric mean
#'
#' @param data an array of numeric values
#' @return \code{gm}
#' @examples
#' data <-  rnorm(1000)
#' geo_mean(data)
#'
#' @import stats utils
#' @export


geo_mean <- function(data) {
    log_data <- log(data)
    gm <- exp(mean(log_data[is.finite(log_data)]))
    return(gm)
}

