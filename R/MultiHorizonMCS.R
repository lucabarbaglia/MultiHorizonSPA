#' @title Multihorizon model confidence set and p-values
#'
#' @description Produces Multi-Horizon Model Confidence Set and p-values as described in section 2.2. of Quadvlieg (2021).
#' @param Losses A list containing M matrices of dimension T by H. Each matrix contains the losses for a forecasting method. Rows correspond to time periods (T rows), columns correspond to forecast horizons (H columns).
#' @param alpha_t Alpha level of critical values for pairwise tests. See Quadvlieg (2021).
#' @param alpha_mcs Alpha level of the critical value for the MCS test. Denoted by alpha tilde in Quadvlieg (2021) section 2.2
#' @param weights the 1 x H vector of weights for the losses at different horizons. For instance \code{weights <- matlab::ones(1,20)/20}
#' @param L integer, the parameter for the moving block bootstrap
#' @param B integer, the number of bootstrap iterations. Default 999
#' @param unif_or_average If equals "a", then the confidence set is calculated based on average SPA. If equals "u", then the confidence set is calculated based on uniform SPA.
#' @export
#' @return Returns a list of length 2. The elements of the list are:
#' \item{MCS_set}{A M by 2 matrix. The first column contains the indices of the methods for which the p-value is greater than or equal to 1-alpha_mcs. The second column contains the p-values that are greater than or equal to 1-alpha_mcs.}
#' \item{p_values}{A vector containing the p-values for all methods.}
#' @examples
#' 
#' 
#' Trow <- 20 
#' H <- 12
#' Mmethods <- 9
#' weights <- rep(1/H,H)
#' 
#' loss_list <- vector(mode = "list", length = Mmethods)
#' 
#' loss_list[[1]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
#' loss_list[[2]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
#' loss_list[[3]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
#' loss_list[[4]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
#' loss_list[[5]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
#' loss_list[[6]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
#' loss_list[[7]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
#' loss_list[[8]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
#' loss_list[[9]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
#' loss_list[[10]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
#' 
#' 
#' MultiHorizonMCS(loss_list, L=3,B=5,unif_or_average = 'u')
#'
#'
#' @export
MultiHorizonMCS <- function(Losses,
                            alpha_t = 0.05,
                            alpha_mcs = 0.05,
                            weights = NULL, 
                            L, 
                            B=999,
                            unif_or_average = "u"){
  
  
  if(!(unif_or_average %in% c("a","u"))){
    stop("unif_or_average must equal 'u' for uniform or 'a' for average.")
  }
  
  
  Set <- Losses
  p_values <- rep(0, length(Losses))
  IDs <- 1:length(Losses)
  
  Trows <- nrow(Losses[[1]])
  Hcols <- ncol(Losses[[1]])
  
  Mmethods <- length(Losses)
  
  if(unif_or_average =="a" &  is.null(weights)){
    weights <- rep(1/Hcols, Hcols)
  }
  
  # print("weights = ")
  # print(weights)
  
  t_uSPA <- matrix(0, nrow = Mmethods, ncol = Mmethods)
  c_uSPA <- matrix(0, nrow = Mmethods, ncol = Mmethods)
  
  for(i in 1:Mmethods){
    for(j in 1:Mmethods){
      if(i!=j){
        LossDiff <- Set[[i]]-Set[[j]]
        
        if(unif_or_average == "u"){
          
          # print('ncol(LossDiff) = ')
          # print(ncol(LossDiff))
          # print('nrow(LossDiff) = ')
          # print(nrow(LossDiff))
          
          bootout <- Bootstrap_uSPA(LossDiff, 
                                    L, 
                                    B)
          
        }else{
          bootout <- Bootstrap_aSPA(LossDiff, 
                                    weights, 
                                    L, 
                                    B)
          
        }
        
        if(unif_or_average == "u"){
          t_uSPA[i,j] <- bootout$t_uSPA
        }else{
          t_uSPA[i,j] <- bootout$t_aSPA
        }
        c_uSPA[i,j] <- quantile(bootout$t_uSPA_b, alpha_t)
      }
    }
  }
  
  
  #create 3 dimensional arrays that contain all the bootstrapped values values
  
  t_uSPA_b <- array(0, c(Mmethods, Mmethods, B));  
  c_uSPA_b <- array(0, c(Mmethods, Mmethods, B));  
  
  for(b in 1:B){
    id <- Get_MBB_ID(Trows, L)
    for(i in 1:Mmethods){
      
      for(j in 1:Mmethods){
        if(i!=j){
          
          if(unif_or_average == "a"){
            # if(is.null(weights)){
            #   Loss_Diff <- Set[[i]]-Set[[j]]
            # }else{
            Loss_Diff <- Set[[i]]-Set[[j]]
            Loss_Diff <- as.matrix(apply(matlab::repmat((weights), T, 1)*Loss_Diff, 1, sum))
            # }
            
          }else{
            Loss_Diff <- Set[[i]]-Set[[j]]
            
          }
          
          
          if(unif_or_average == "a"){
            Demeaned_Loss_Diff <- (Loss_Diff - matlab::repmat(mean(Loss_Diff),Trows,1))[id,]
          }else{
            Demeaned_Loss_Diff <- sweep(Loss_Diff,2,
                                        apply(Loss_Diff, 2, mean))[id,]
          }
          
          
          
          # Weighted_Loss_Diff <- as.matrix(apply(matlab::repmat(t(weights), Trows, 1)*Loss_Diff, 1, sum))
          # d_ij <- mean(Weighted_Loss_Diff)
          if(unif_or_average == "a"){
            d_ij <- mean(Demeaned_Loss_Diff)
            t_uSPA_temp <- as.numeric(sqrt(Trows)*d_ij/sqrt(QS(Demeaned_Loss_Diff)))
            
          }else{
            d_ij <- apply(Demeaned_Loss_Diff, 2, mean)
            t_uSPA_temp <- min(sqrt(Trows)*d_ij/sqrt(QS(Demeaned_Loss_Diff))) 
          }
  
          
          #applying bootstrap again
          t_uSPA_b_temp <- matlab::zeros(B, 1)
          # Demeaned_Loss_Diff <- Loss_Diff- matlab::repmat(d_ij, Trows, 1)
          
          for (b2 in 1:B){
            id2 <- Get_MBB_ID(Trows, L)
            b_lossdiff <- as.matrix(Demeaned_Loss_Diff[id2, ])
            omega_b <- MBB_Variance(b_lossdiff,L)
            
            if(unif_or_average == "a"){
              t_uSPA_b_temp[b2] <- as.numeric(sqrt(Trows)*mean(b_lossdiff)/sqrt(omega_b))
            }else{
              t_uSPA_b_temp[b2] <- min(sqrt(Trows)*apply(b_lossdiff,2,mean)/sqrt(omega_b))
            }
          }
          
          
          
          t_uSPA_b[i,j,b] <- t_uSPA_temp
          c_uSPA_b[i,j,b] <- quantile(t_uSPA_b_temp, alpha_t)
          
        }# end if i!=j
      } #end loop over j
    }#end loop over i
  }#end loop over b
  
  
  #now implement while loop over M
  
  while(Mmethods > 0){
    
    
    # print("Mmethods=")
    # print(Mmethods)
    # 
    # print("t_uSPA_b=")
    # print(t_uSPA_b)
    
    t_max_uSPA 	<- max(t_uSPA	- c_uSPA)
    if(Mmethods==1){
      t_max_uSPA_b <- max(t_uSPA_b- c_uSPA_b, na.rm = TRUE)
    }else{
      t_max_uSPA_b <- apply(t_uSPA_b - c_uSPA_b,MARGIN = c(3), FUN = max, na.rm = TRUE)
      
    }

    # print("t_max_uSPA_b=")
    # print(t_max_uSPA_b)
    

    crits <- t_uSPA	- c_uSPA - 9999999*diag(Mmethods)
    index <- which.max(apply(crits,1,max, na.rm = TRUE))
    if(length(index)==0){index <- 1} #replicating code. Perhaps should throw error here?
    
    
    # print("Mmethods=")
    # print(Mmethods)
    
    # if(Mmethods==1){print("Mmethods==1")}
    # if(ncol(t_max_uSPA_b)==1){print("ncol(t_max_uSPA_b)==1")}
    # if(is.vector(t_max_uSPA_b)){print("t_max_uSPA_b is a vector")}
    
    
    # print("t_max_uSPA_b=")
    # print(t_max_uSPA_b)
    # 
    # print("t_max_uSPA=")
    # print(t_max_uSPA)
    
    
    # if(Mmethods==1){
      p_temp <- mean((t_max_uSPA < t_max_uSPA_b))
    # }else{
    #   p_temp <- apply((t_max_uSPA < t_max_uSPA_b),2,mean, na.m = TRUE)
    #   
    # }
    
    # print("p_temp=")
    # print(p_temp)
    # 
    # print("p_values=")
    # print(p_values)
    # 
    # 
    # print("index=")
    # print(index)
    
    
    #CHECK THIS LINE
    p_values[IDs[index]] <- max(p_temp, p_values,na.rm = TRUE)   #apply(p_temp,max(p_values),1,max);
    
    #if the following, then would not give one element
    #but then why are row maximums used in the original code?
    #p_values[IDs[index]] = apply(cbind( p_temp,rep(max(p_values), length(p_temp))  ),1,max);
    
    if(length(IDs)==1){
      p_values[IDs[index]] <- 1
    }
    
    # models <- array(NA, c(Mmethods, Mmethods, 2))
    # 
    # for(i in 1:Mmethods){
    #   for(j in 1:Mmethods){
    #     models[i,j,1] <- i
    #     models[i,j,2] <- j
    #   }
    # }
    
    IDs <- IDs[ - index ]
    #this line is unnecessary?
    # Set <- Set[-index,]
    
    
    
    # deletes <- (models[,,1] == index) + (models[,,2] == index);
    
    if(Mmethods == 1){
      
    }else{
      t_uSPA <- t_uSPA[-index, - index]
      c_uSPA <- c_uSPA[-index, - index]
      t_uSPA_b <- t_uSPA_b[-index, - index,]
      c_uSPA_b <- c_uSPA_b[-index, - index,]
      
    }
    
    
    
    Mmethods <- Mmethods -1
    
    
    
  }
  
  
  MCS_set <- cbind(1:length(Losses), p_values)
  
  MCS_set <- MCS_set[ !( p_values < 1-alpha_mcs)  ,]
  
  
  ret_list <- list()
  
  ret_list$MCS_set <- MCS_set
  ret_list$p_values <- p_values
  
  return(ret_list)
  
}

