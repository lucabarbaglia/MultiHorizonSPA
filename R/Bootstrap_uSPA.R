#' Bootstrap uSPA
#' 
#' Implements Bootstrap Algorithm 1 for the test for Uniform SPA of Quaedvlieg (2021).
#' Bootstrap is based on a moving block bootstrap with length L.
#'  
#' @keywords internal
#' @noRd

Bootstrap_uSPA <- function(Loss_Diff, L, B){
  # Implements Bootstrap Algorithm 1 for the test for Uniform SPA of Quaedvlieg (2018).
  # Bootstrap is based on a moving block bootstrap with length L.
  
  #TT <- matlab::size(Loss_Diff, 1)
  TT <- nrow(Loss_Diff)
  
  d_ij <- apply(Loss_Diff, 2, mean)
  t_uSPA <- min(sqrt(TT)*d_ij/sqrt(QS(Loss_Diff))) 
  
  t_uSPA_b <- matlab::zeros(B, 1)
  Demeaned_Loss_Diff <- Loss_Diff- matlab::repmat(d_ij, TT, 1)
  
  
  
  for (b in 1:B){
    id <- Get_MBB_ID(TT, L)
    b_lossdiff <- Demeaned_Loss_Diff[id, ]
    omega_b <- MBB_Variance(b_lossdiff,L)
    t_uSPA_b[b] <- min(sqrt(TT)*apply(b_lossdiff,2,mean)/sqrt(omega_b))
  }
  return(list("t_uSPA" = t_uSPA, "t_uSPA_b" = t_uSPA_b))
}

