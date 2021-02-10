#' QS
#' 
#' Returns column-wise QS HAC estimator of the variance
#'  
#' @keywords internal
#' @noRd

QS <- function(y){
    TT <- nrow(y)
    N <- ncol(y)
    bw <- 1.3*TT^(1/5)
    weight <- QS_Weights((1:(TT-1))/bw)
    omega <- matlab::zeros(1,N)
    for (i in 1:N){
        workdata <- y[,i]-mean(y[,i])
        omega[i] <-  t(workdata)%*%workdata/TT
        for (j in 1:(TT-1)){
            omega[i] <- omega[i] + 2*weight[j]*(workdata[1:(TT-j)]%*%workdata[(j+1):TT])/TT
        }
    }
    return(omega)
}

