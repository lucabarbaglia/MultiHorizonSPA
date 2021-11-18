#' Fast Test average Superior Predictive Ability
#' 
#' Implements the test for average Superior Predictive Ability (aSPA) of Quaedvlieg (2021)
#' @param LossDiff the T x H matrix forecast path loss differential
#' @param weights the 1 x H vector of weights for the losses at different horizons. For instance \code{weights <- matlab::ones(1,20)/20}
#' @param L integer, the parameter for the moving block bootstrap
#' @param B integer, the number of bootstrap iterations. Default 999
#' @return A list containing two objects:
#' \item{"p_value"}{the p-value for aSPA}
#' \item{"t_aSPA"}{the statistics for aSPA}
#' @references Quaedvlieg, Rogier. "Multi-horizon forecast comparison." Journal of Business & Economic Statistics 39.1 (2021): 40-53.
#' @seealso \code{\link{Test_uSPA}}
#' @examples
#' ## Test for aSPA and uSPA
#' Trow <- 200 
#' H <- 12
#' Mmethods <- 5
#' weights <- rep(1/H,H)
#' Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)
#' 
#' Fast_Test_aSPA(LossDiff=Losses, weights=weights, L=3, B=10)
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export

Fast_Test_aSPA <- function(LossDiff, 
                           weights = NULL, 
                           L, 
                           B=999){
  
  
  if (!is.matrix(LossDiff)){LossDiff <- as.matrix(LossDiff)}
  
  Hcols <- ncol(LossDiff)
  
  Mmethods <- length(Losses)
  
  if(is.null(weights)){
    weights <- rep(1/Hcols, Hcols)
  }
  if(length(weights)!=Hcols){
    print("length of weights not equal to number of columns (horizons), setting equal weights")
    weights <- rep(1/Hcols, Hcols)
  }
  
  
  ret <- Test_aSPA_cpp(as.matrix(LossDiff),
                       weights,
                       L, 
                       B)
  
  # bootout <- Bootstrap_uSPA(LossDiff, L, B)
  # p_value <- mean(bootout$t_uSPA < bootout$t_uSPA_b)
  # 
  
  
  
  return(list("p_value"=ret[[0]], "t_aSPA"=ret[[1]]))
  
  
}

