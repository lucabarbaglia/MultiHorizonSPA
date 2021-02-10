#' Get_MBB_ID
#' 
#' Get_MBB_ID obtains ids of re-sampled observations using a moving block bootstrap with blocks of length L
#'  
#' @keywords internal
#' @noRd


Get_MBB_ID <- function(TT, L){
  #Get_MBB_ID obtains ids of resampled observations using a
  #moving block bootstrap with blocks of length L
  id <- matlab::zeros(TT,1)
  id[1] <- matlab::ceil(TT*stats::runif(1))
  for (t in 2:TT){
    if (matlab::rem(t,L) == 0){
      id[t] <- matlab::ceil(TT*stats::runif(1))
    } else {
      id[t] <- id[t-1]+1
    }
    if (id[t] > TT){
      id[t] <- 1
    }
  }
  return(id)
}

