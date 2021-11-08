#' Bootstrap aSPA
#' 
#' Implements Bootstrap Algorithm 1 for the test for Average SPA of Quaedvlieg (2021).
#' Bootstrap is based on a Moving Block Bootstrap with parameter L.
#' Weights are the weights on the various horizons, size conforming to Loss_Diff.
#'  
#' @keywords internal
#' @noRd

Bootstrap_aSPA <- function(Loss_Diff, weights, L, B){
  TT <- nrow(Loss_Diff)
  
  Weighted_Loss_Diff <- as.matrix(apply(matlab::repmat((weights), TT, 1)*Loss_Diff, 1, sum))
  d_ij <- mean(Weighted_Loss_Diff)
  t_aSPA <- as.numeric(sqrt(TT)*d_ij/sqrt(QS(Weighted_Loss_Diff)))
  
  t_aSPA_b <- matlab::zeros(B,1)
  Demeaned_Loss_Diff <- Weighted_Loss_Diff - matlab::repmat(d_ij,TT,1)
  for (b in 1:B){
    id <- Get_MBB_ID(TT, L)
    b_lossdiff <- as.matrix(Demeaned_Loss_Diff[id,])
    zeta_b <- MBB_Variance(b_lossdiff,L)
    t_aSPA_b[b] <- sqrt(TT)*( mean(b_lossdiff))/sqrt(zeta_b)
  }
  return(list("t_aSPA" = t_aSPA, "t_aSPA_b" = t_aSPA_b))
}






