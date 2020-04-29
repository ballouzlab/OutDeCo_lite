#' Recurrence count
#'
#' @param data A matrix
#' @param nmax number of datasets
#' @return \code{res}
#' @examples
#'
#' data <-  rowSums(matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100))
#' count_recur(as.matrix(data), 10)
#'
#' @import plyr stats utils
#' @export
#'
#'

count_recur <- function(data, nmax) {
    freq <- plyr::count(data[, 1])
    res <- matrix(0, nrow = nmax, ncol = 1)
    rownames(res) <- 0:(nmax - 1)
    m <- match(0:(nmax - 1), freq[, 1])
    f.r <- !is.na(m)
    f.f <- m[f.r]
    res[f.r, 1] <- freq[f.f, 2]
    return(res)
}
