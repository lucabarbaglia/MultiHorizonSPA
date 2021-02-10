#' Test uniform Superior Predictive Ability
#' 
#' Implements the test for uniform Superior Predictive Ability (uSPA) of Quaedvlieg (2021)
#' @param LossDiff the TxH matrix forecast path loss differential
#' @param weights the 1xH vector of weights for the losses at different horizons. For instance \code{weights <- matlab::ones(1,20)/20}
#' @param L the parameter for the moving block bootstrap
#' @return A list containg two objects:
#' \item{"p_value"}{the p-value for uSPA}
#' \item{"t_uSPA"}{the statistics for uSPA}
#' @references Quaedvlieg, Rogier. "Multi-horizon forecast comparison." Journal of Business & Economic Statistics 39.1 (2021): 40-53.
#' @seealso \code{\link{Test_aSPA}}
#' @examples
#' ## Test for uSPA
#' data(LossDiff_uSPA)
#' Test_uSPA(LossDiff=LossDiff_uSPA, L=3)
#' 
#' @author Luca Barbaglia \email{luca.barbaglia@ec.europa.eu}
#' 
#' @export

Test_uSPA <- function(LossDiff, L){
  if (!is.matrix(LossDiff)){LossDiff <- as.matrix(LossDiff)}
  bootout <- Bootstrap_uSPA(LossDiff, L)
  p_value <- mean(bootout$t_uSPA < bootout$t_uSPA_b)
  return(list("p_value"=p_value, "t_uSPA"=bootout$t_uSPA))
}

