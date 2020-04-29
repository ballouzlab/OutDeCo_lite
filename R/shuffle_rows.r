#' Shuffle rows
#'
#' @param data A matrix
#' @return \code{output}
#' @examples
#'
#' data <- matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' shuffle_rows(data)
#'
#'
#' @export
#'

shuffle_rows <- function(data) {
    nc <- dim(data)[2]
    nr <- dim(data)[1]
    output <- sapply(1:nr, function(i) data[i, sample(nc)])
    return(output)
}
