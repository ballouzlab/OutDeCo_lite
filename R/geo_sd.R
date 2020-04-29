#' geometric standard deviation
#'
#' @param data an array of numeric values
#' @return \code{gs}
#' @examples
#' data <-  rnorm(1000)
#' geo_sd(data)
#'
#' @import stats utils
#' @export
#'


geo_sd <- function(data) {
    log_data <- log(data)
    gs <- exp(sd(log_data[is.finite(log_data)]))
    return(gs)
}
