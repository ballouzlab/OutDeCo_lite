#' Probability calc
#'
#' @param data A matrix
#' @param pp significance/FDR threshold
#' @param nr number of repeats to run
#' @return \code{res}
#' @examples
#'
#' genesets <- matrix(  rbinom(1000,1,0.05), ncol=10, nrow=100)
#' fdrs_recur <- calc_fdrs_recur(genesets)
#'
#' @import stats utils
#' @export
#'


calc_fdrs_recur <- function(data, pp = 0.05, nr = 1000) {
    
    temp1 <- lapply(1:nr, function(i) shuffle_cols(data))
    temp2 <- sapply(1:nr, function(i) rowSums(temp1[[i]], na.rm = T))
    
    nmax <- dim(data)[2] + 1
    
    recur <- rowSums(data, na.rm = T)
    ob <- count_recur(as.matrix(recur), nmax)
    ex <- rowSums(sapply(1:nr, function(j) count_recur(as.matrix(temp2[, j]), nmax)))/nr
    
    test <- cbind(ob, ex)
    
    observed <- cbind((rev(cumsum(rev(test[, 1])))), 0:(nmax - 1))
    expected <- cbind((rev(cumsum(rev(test[, 2])))), 0:(nmax - 1))
    
    FDR <- expected[, 1]/observed[, 1]
    Pt <- expected[min(which(FDR < pp)), 2]
    sig <- sum(test[(which(FDR < pp)), 1], na.rm = T)
    res <- list(FDR, test, Pt, sig, pp)
    names(res) <- c("FDR", "test", "Pt", "sig", "pp")
    return(res)
}

