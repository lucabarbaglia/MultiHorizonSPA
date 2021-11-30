#' Fast Test uniform Superior Predictive Ability
#' 
#' Implements the test for uniform Superior Predictive Ability (uSPA) of Quaedvlieg (2021)
#' @param LossDiff the T x H matrix forecast path loss differential
#' @param L the parameter for the moving block bootstrap
#' @param B integer, the number of bootstrap iterations. Default 999
#' @return A list containing two objects:
#' \item{"p_value"}{the p-value for uSPA}
#' \item{"t_uSPA"}{the statistics for uSPA}
#' @references Quaedvlieg, Rogier. "Multi-horizon forecast comparison." Journal of Business & Economic Statistics 39.1 (2021): 40-53.
#' @seealso \code{\link{Test_aSPA}}
#' @examples
#' ## Test for uSPA
#' Trow <- 200 
#' H <- 12
#' Mmethods <- 5
#' Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)
#' 
#' Fast_Test_uSPA(LossDiff=Losses, L=3, B=10)
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export

Fast_Test_uSPA <- function(LossDiff, 
                           L, 
                           B=999,
                           num_cores = 1,
                           seed = runif(1, 0, .Machine$integer.max)){
  
  
  if (!is.matrix(LossDiff)){LossDiff <- as.matrix(LossDiff)}
  
  ret <- Test_uSPA_cpp(as.matrix(LossDiff), 
                       L, 
                       B,
                       num_cores,
                       seed)
  
  # bootout <- Bootstrap_uSPA(LossDiff, L, B)
  # p_value <- mean(bootout$t_uSPA < bootout$t_uSPA_b)
  # 
  
  
  
  return(list("p_value"=ret[[1]], "t_uSPA"=ret[[2]]))
  
  
}

