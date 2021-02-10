#' Test average Superior Predictive Ability
#' 
#' Implements the test for average Superior Predictive Ability (aSPA) of Quaedvlieg (2021)
#' @param LossDiff the TxH matrix forecast path loss differential
#' @param weights the 1xH vector of weights for the losses at different horizons. For instance \code{weights <- matlab::ones(1,20)/20}
#' @param L the parameter for the moving block bootstrap
#' @return A list containg two objects:
#' \item{"p_value"}{the p-value for aSPA}
#' \item{"t_aSPA"}{the statistics for aSPA}
#' @references Quaedvlieg, Rogier. "Multi-horizon forecast comparison." Journal of Business & Economic Statistics 39.1 (2021): 40-53.
#' @seealso \code{\link{Test_uSPA}}
#' @examples
#' ## Test for aSPA and uSPA
#' data(LossDiff_aSPA)
#' weights <- matlab::ones(1,20)/20
#' Test_aSPA(LossDiff=LossDiff_aSPA, weights=weights, L=3)
#' Test_uSPA(LossDiff=LossDiff_aSPA, L=3)
#' 
#' @author Luca Barbaglia \email{luca.barbaglia@ec.europa.eu}
#' 
#' @export



Test_aSPA <- function(LossDiff, weights, L){
  if (!is.matrix(LossDiff)){LossDiff <- as.matrix(LossDiff)}
  bootout <- Bootstrap_aSPA( LossDiff, weights, L)
  p_value <- mean(bootout$t_aSPA < bootout$t_aSPA_b)
  return(list("p_value"=p_value, "t_aSPA"=bootout$t_aSPA))
  
}


