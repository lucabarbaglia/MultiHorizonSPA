#' Test uniform Superior Predictive Ability
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
#' data(LossDiff_uSPA)
#' Test_uSPA(LossDiff=LossDiff_uSPA, L=3, B=10)
#' 
#' @author Luca Barbaglia \url{https://lucabarbaglia.github.io/}
#' 
#' @export

Test_uSPA <- function(LossDiff, L, B=999){
  if (!is.matrix(LossDiff)){LossDiff <- as.matrix(LossDiff)}
  bootout <- Bootstrap_uSPA(LossDiff, L, B)
  p_value <- mean(bootout$t_uSPA < bootout$t_uSPA_b)
  return(list("p_value"=p_value, "t_uSPA"=bootout$t_uSPA))
}

