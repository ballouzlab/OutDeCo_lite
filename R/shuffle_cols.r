#' Shuffle columns
#'
#' @param data A matrix
#' @return \code{output}
#' @examples
#' data <- matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' shuffle_cols(data)
#'
#' @export
#'

shuffle_cols <- function(data) {
    nc <- dim(data)[2]
    nr <- dim(data)[1]
    output <- (sapply(1:nc, function(i) data[sample(nr), i]))
    return(output)
}

