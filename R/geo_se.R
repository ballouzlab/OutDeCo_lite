#' geometric standard error
#'
#' @param data an array of numeric values
#' @return \code{gse}
#' @examples
#'
#' data <-  rnorm(1000)
#' geo_se(data)
#'
#' @import stats utils
#' @export

geo_se <- function(data) {
    gs <- geo_sd(data)
    log_data <- log(data)
    gse <- gs/sqrt(sum(is.finite(log_data)))
    return(gse)
}

