#' Recurrence to frequency matrix
#'
#' @param recurs A vector of gene recurrences
#' @return \code{mat10}
#' @examples
#'
#' data_pre <- as.matrix(rowSums(matrix(  rbinom(1000,1,0.1), ncol=10, nrow=100)))
#' data_post <- data_pre - (data_pre * rbinom(100, 1, 0.2 ))
#' recurs = cbind(data_post, data_pre)
#'
#' get_recur_mat(recurs)
#'
#' @import plyr stats utils graphics
#' @export
#'

get_recur_mat <- function(recurs) {
    temp1 <- recurs
    temp <- plyr::count(temp1)
    x <- "x"
    y <- "y"
    freq <- "freq"
    colnames(temp) <- c(x, y, freq)
    temp.mat <- tidyr::spread(temp, key = x, value = freq)
    rownames(temp.mat) <- temp.mat[, 1]
    temp.mat <- temp.mat[, -1]
    temp.mat[is.na(temp.mat)] <- 0
    temp.mat <- as.matrix(temp.mat)
    
    mat10 <- log10(temp.mat) + 1
    mat10[!is.finite(mat10)] <- 0
    return(mat10)
}
