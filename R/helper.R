#' OutDeCo helper functions
#'

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



