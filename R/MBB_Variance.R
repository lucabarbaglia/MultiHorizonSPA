#' MBB_Variance
#' 
#' Computes the 'natural variance' estimator of the re-sampled data
#'  
#' @keywords internal
#' @noRd


MBB_Variance <- function(y, L){
	TT <- nrow(y)
	N <- ncol(y)
	omega <- matlab::zeros(1,N)
	ydem <- y - matlab::repmat(apply(y, 2, mean),TT,1)
    K <- floor(TT/L)
    for (n in 1:N){
        temp <- matlab::reshape(as.matrix(ydem[1:(K*L),n]), c(K,L))
        omega[n] <- mean(apply(temp, 2, sum)^2)/L
    }
    return(omega)
}
