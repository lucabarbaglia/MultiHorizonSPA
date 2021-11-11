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
    
    # print("bw = ")
    # print(bw)
    
    weight <- QS_Weights((1:(TT-1))/bw)
    omega <- matlab::zeros(1,N)
    
    # print("weight = ")
    # print(weight)
    # 
    
    for (i in 1:N){
        workdata <- y[,i]-mean(y[,i])
        omega[i] <-  t(workdata)%*%workdata/TT
        
        # print("workdata = ")
        # print(workdata)
        # 
        # print("omega[i] = ")
        # print(omega[i])
        
        for (j in 1:(TT-1)){
            
            # print("workdata[1:(TT-j)] = ")
            # print(workdata[1:(TT-j)])
            
            
            omega[i] <- omega[i] + 2*weight[j]*(workdata[1:(TT-j)]%*%workdata[(j+1):TT])/TT
        }
    }
    
    # print("omega = ")
    # print(omega)
    
    return(omega)
}

