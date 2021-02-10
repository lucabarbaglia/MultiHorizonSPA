#' QS_Weights
#' 
#' QS_WEIGHTS, from HAC() function
#'  
#' @keywords internal
#' @noRd

QS_Weights <- function(x){
    argQS <- 6*pi*x/5
    w1 <- 3/(argQS^2)
    w2 <- (sin(argQS)/argQS)-cos(argQS)
    wQS <- w1*w2
    wQS[wQS == 0] <- 1
    return(wQS)
}

