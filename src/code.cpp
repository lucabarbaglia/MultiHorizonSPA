



# include <RcppArmadillo.h>
#include <cmath>

// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec QS_Weights_cpp(arma::vec x){
  
  arma::vec argQS = 6*arma::datum::pi*x/5;
  arma::vec w1 = 3/( arma::pow(argQS ,2));  //note: arma::pow() should make element-wise calculations
  arma::vec w2 = (arma::sin(argQS)/argQS)-arma::cos(argQS);
  arma::vec wQS = w1%w2;
  wQS.replace(0,1);
  //alternative to above line
  //wQS.elem( find(wQS == 0) ).ones();
  
  return(wQS);
}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec QS_cpp(arma::mat y){
  
  
  int TT = y.n_rows;
  int N  = y.n_cols;
  
  double bw = 1.3*std::pow(double(TT), double(0.2));
  
  // Rcout << "bw = " << bw << ". \n";
  
  
  
  arma::vec weightarg = arma::regspace<arma::vec>(1,TT-1)/bw;
  arma::vec weight = QS_Weights_cpp( weightarg );
  arma::vec omega = arma::zeros<arma::vec>(N);
    
  // Rcout << "Line 50 . \n";
    
  
  for(int i=0; i<N;i++){
    arma::vec workdata = y.col(i) - arma::mean(y.col(i) )  ;
    double temp = arma::as_scalar(workdata.t()*workdata) ;

    
    
    omega(i) = temp / TT ;
    
    // Rcout << " workdata = " <<  workdata << ". \n";
    // 
    // Rcout << " omega(i)  = " <<  omega(i)  << ". \n";
    
    
    for(int j=0; j<TT-1 ;j++){
      arma::uvec inds1 = arma::regspace<arma::uvec>(0,TT-j-2); 
      arma::uvec inds2 = arma::regspace<arma::uvec>(j+1,TT-1);
      
      arma::vec tempvec1 = workdata(inds1);
      arma::vec tempvec2 = workdata(inds2);
        
        
        // Rcout << " tempvec1  = " <<  tempvec1  << ". \n";
        
      omega(i) +=  2*(weight(j))*arma::as_scalar(tempvec1.t()*tempvec2)/double(TT)  ;
    }
  }
  
  // Rcout << "omega = " << omega << ". \n";
  // Rcout << "Line 83 . \n";
  
  
  return(omega);
}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]

#include <xoshiro.h>
// #include <dqrng_distribution.h>
// #include <dqrng.h>

// [[Rcpp::export]]
arma::uvec Get_MBB_ID_cpp(int TT, int L){
  
  arma::uvec id = arma::zeros<arma::uvec>(TT);
  
  
  std::random_device device;

  dqrng::xoshiro256plus gen(device());              // properly seeded rng
  // dqrng::xoshiro256plus gen(seed);              // properly seeded rng
  
  
  
  std::uniform_real_distribution<> dis_cont_unif(0, 1);
  
  double temprand = dis_cont_unif(gen);
  
  id(0) = std::ceil((double)TT*temprand)-1;

  for(int t=1; t<TT;t++){
    if( (t+1) % L == 0){
      id(t) = std::ceil(TT*dis_cont_unif(gen))-1;
    }else{
      id(t) = id(t-1)+1;
    }
    
    if(id(t) > TT-1){
      id(t) = 0 ;
    }
  
  }//end of for loop over t
  
  return(id);
}


//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec MBB_Variance_cpp(arma::mat y, int L){
  
  double TT = y.n_rows;
  int N  = y.n_cols;
  
  arma::vec omega = arma::zeros<arma::vec>(N);
  
  
  // arma::vec colmeans(N) ;
  // 
  // for(int i=0; i<N;i++){
  //   colmeans(i) = arma::mean(y.col(i));
  // }
  // arma::mat ydem = y.each_row() - colmeans  ; 
  
  arma::mat ydem = y.each_row() - arma::mean(y,0)  ; 
  
  int K = std::floor(TT/L);
  
  arma::uvec inds1 = arma::regspace<arma::uvec>(0,K*L-1); 
  
  
  for(int n=0; n<N;n++){
    
    arma::vec coltemp1 = ydem.col(n);
    arma::vec coltemp2 = coltemp1.elem(inds1);
    
    arma::mat temp = arma::reshape(coltemp2, K, L);
    
    // arma::vec temp2 = arma::sum(temp)

    omega(n) = arma::mean(arma::pow(arma::sum(temp),2))/double(L);
      
  }
  
  
  return(omega);
}



//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]

#include <xoshiro.h>
// #include <dqrng_distribution.h>
// #include <dqrng.h>

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

// [[Rcpp::export]]
arma::field<arma::vec> Bootstrap_aSPA_cpp(arma::mat Loss_Diff,arma::vec weights, int  L, int B,
                                          int ncores,
                                          int seed){

    
  double TT = Loss_Diff.n_rows;
  
  // % is element-wise multiplication
  //row sums after multiplying each row element by weights (horizon weights)
  
  //really a vector, saving as matrix for similarity to R code
  arma::mat Weighted_Loss_Diff = arma::sum((Loss_Diff.each_row() % weights.t()),1);
  
  //instead of column means, appears to be one mean of weighted row sums
  double d_ij = arma::as_scalar(mean(Weighted_Loss_Diff));
  
  
  
  //vector of length one. Keeping as vector for field output
  arma::vec t_aSPA = std::sqrt(TT)*d_ij/arma::sqrt(QS_cpp(Weighted_Loss_Diff));
  
  //double t_aSPA = arma::as_scalar(std::sqrt(TT)*d_ij/arma::sqrt(QS_cpp(Weighted_Loss_Diff)));
  
  arma::vec t_aSPA_b = arma::zeros<arma::vec>(B);
  
  // demean the vector (saved as matrix)
  arma::mat Demeaned_Loss_Diff = Weighted_Loss_Diff - d_ij;
  
  std::random_device device;
  
  // dqrng::xoshiro256plus gen(device());              // properly seeded rng
  dqrng::xoshiro256plus gen(seed); 
  
  std::uniform_real_distribution<> dis_cont_unif(0, 1);
  
  //eventually parallelize this using number of cores equal to the minimum of the number available and B
  int ncoretemp = std::min(ncores, B);
  
  // Rcout << "line 758. \n" ;
  
  
#pragma omp parallel num_threads(ncoretemp)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps
  // 
#pragma omp for
  for(int b=0; b<B;b++){
    //arma::uvec id = Get_MBB_ID_cpp(TT,L );
    
    arma::uvec id = arma::zeros<arma::uvec>(TT);
    
    
    // std::random_device device;
    // 
    // dqrng::xoshiro256plus gen(device());              // properly seeded rng
    //dqrng::xoshiro256plus gen(seed);              // properly seeded rng
    
    
    
    
    double temprand = dis_cont_unif(lgen);
    
    id(0) = std::ceil((double)TT*temprand)-1;
    
    for(int t=1; t<TT;t++){
      if( (t+1) % L == 0){
        id(t) = std::ceil(TT*dis_cont_unif(lgen))-1;
      }else{
        id(t) = id(t-1)+1;
      }
      
      if(id(t) > TT-1){
        id(t) = 0 ;
      }
      
    }//end of for loop over t
    
    
    
    
    
    // Rcout << "id = " << id << ". \n";
    
    arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id);
    
    // Rcout << "line 260 . \n";
    
    arma::vec zeta_b = MBB_Variance_cpp(b_lossdiff, L);
    //only one column, so arma::mean gives vector of length one
    t_aSPA_b(b) = arma::as_scalar(std::sqrt(TT)*( arma::mean(b_lossdiff))/arma::sqrt(zeta_b));
  }
 
}//end pragma code   
  
  
  arma::field<arma::vec> ret_f(2);
    
  ret_f(0) = t_aSPA;
  ret_f(1) = t_aSPA_b;
  
  return(ret_f);
}



//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]

#include <xoshiro.h>
// #include <dqrng_distribution.h>
// #include <dqrng.h>

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

// [[Rcpp::export]]
arma::field<arma::vec> Bootstrap_uSPA_cpp(arma::mat Loss_Diff, int  L, int B,
                                          int ncores,
                                          int seed){
  
  // Rcout << "line 290 . \n";
  
  double TT = Loss_Diff.n_rows;
  
  // Rcout << "line 294 . \n";
  
  //vector of column means
  arma::rowvec d_ij = arma::mean(Loss_Diff);
  
  // Rcout << "line 296 . \n";
  
  //should output a vector of length 1
  //numerator is vector of length H (num columns)
  //denominator is vector of length H (num columns)
  // element-wise division, then minimum over vector
  
  //perhaps arma::min returns a double
  //arma::vec t_uSPA = arma::min(std::sqrt(TT)*d_ij/arma::sqrt(QS_cpp(Loss_Diff))) ;
  
  arma::vec t_uSPA(1);
  t_uSPA(0) = arma::min(std::sqrt(TT)*d_ij.t()/arma::sqrt(QS_cpp(Loss_Diff))) ;

  arma::vec t_uSPA_b = arma::zeros<arma::vec>(B);
  
  
  // Rcout << "line 315 . \n";
  
  arma::mat Demeaned_Loss_Diff = Loss_Diff.each_row() - d_ij;
  
  // Rcout << "line 319 . \n";
  
  std::random_device device;
  
  // dqrng::xoshiro256plus gen(device()); 
  dqrng::xoshiro256plus gen(seed); 
  
  std::uniform_real_distribution<> dis_cont_unif(0, 1);
  
  int ncoretemp = std::min(ncores, B);
  
  
#pragma omp parallel num_threads(ncoretemp)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps
  // 
#pragma omp for
  for(int b=0; b<B;b++){
    // arma::uvec id = Get_MBB_ID_cpp(TT,L);
    
    arma::uvec id = arma::zeros<arma::uvec>(TT);
    
    
    // std::random_device device;
    // 
    // dqrng::xoshiro256plus gen(device());              // properly seeded rng
    //dqrng::xoshiro256plus gen(seed);              // properly seeded rng
    
    
    double temprand = dis_cont_unif(lgen);
    
    id(0) = std::ceil((double)TT*temprand)-1;
    
    for(int t=1; t<TT;t++){
      if( (t+1) % L == 0){
        id(t) = std::ceil(TT*dis_cont_unif(lgen))-1;
      }else{
        id(t) = id(t-1);
      }
      
      if(id(t) > TT-1){
        id(t) = 0 ;
      }
      
    }//end of for loop over t
    
    // Rcout << "line 350 . \n";
    
    arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id);
    arma::vec omega_b = MBB_Variance_cpp(b_lossdiff, L);
    
    // Rcout << "line 355 . \n";
    
    t_uSPA_b(b) = arma::min(std::sqrt(TT)*( arma::mean(b_lossdiff).t())/arma::sqrt(omega_b));
    
    
  }
  
}//end of pragma code

  // Rcout << "line 359 . \n";
  
  
  arma::field<arma::vec> ret_f(2);
  
  ret_f(0) = t_uSPA;
  ret_f(1) = t_uSPA_b;
  
  return(ret_f);
}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Test_aSPA_cpp(NumericMatrix LossDiff,
                   NumericVector weights,
                   int L,
                   int B,
                   int ncores,
                   int seed){
  

  arma::mat LossDiff_arma = as<arma::mat>(LossDiff);
  arma::vec weights_arma = as<arma::vec>(weights);
  
  arma::field<arma::vec> bootout = Bootstrap_aSPA_cpp(LossDiff_arma, weights_arma, L, B, ncores, seed) ; 
  
  arma::vec temp_vec = bootout(0);
  double t_aSPA = temp_vec(0) ;
  
  arma::vec t_aSPA_b = bootout(1);
  
  arma::uvec temp =  (t_aSPA < t_aSPA_b );
  
  // Rcout << "temp = " << temp << ". \n";
  
  
  //calculating the meaof binary numbers
  double p_value = double(arma::sum(temp))/double(t_aSPA_b.n_elem);
  
  
  
  List ret(2);
  
  ret[0] = p_value;
  ret[1] = t_aSPA;
  
  
  return(ret);
}





//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Test_uSPA_cpp(NumericMatrix LossDiff,
                   int L,
                   int B,
                   int ncores,
                   int seed){
  
  // Rcout << "line 426 . \n";
  
  arma::mat LossDiff_arma = as<arma::mat>(LossDiff);
  
  // Rcout << "line 429 . \n";
  
  arma::field<arma::vec> bootout = Bootstrap_uSPA_cpp(LossDiff_arma, L, B, ncores, seed) ; 
  
  arma::vec temp_vec = bootout(0);
  double t_uSPA = temp_vec(0) ;
  
  arma::vec t_uSPA_b = bootout(1);
  
  arma::uvec temp =  (t_uSPA < t_uSPA_b );
  
  // Rcout << "temp = " << temp << ". \n";
  
  
  //calculating the meaof binary numbers
  double p_value = double(arma::sum(temp))/double(t_uSPA_b.n_elem);
  
  
  
  List ret(2);
  
  ret[0] = p_value;
  ret[1] = t_uSPA;
  
  
  return(ret);
}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]

#include <xoshiro.h>
// #include <dqrng_distribution.h>
// #include <dqrng.h>

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

// [[Rcpp::export]]
List MultiHorizonMCS_cpp(List Losses, 
                         double alpha_t,
                         double alpha_mcs,
                         NumericVector weights,
                         int L, 
                         int B,
                         int unif_or_avg, //1 means uniform, 0 means average
                         int ncores,
                         int seed
                           ){
  

  // Rcout << "line 490. \n" ;
  
  arma::vec alpha_t_vec = { alpha_t };
  arma::vec alpha_mcs_vec = { alpha_mcs };
  
  
  // arma::mat LossDiff_arma = as<arma::mat>(LossDiff);
  arma::vec weights_arma = as<arma::vec>(weights);
  
  // Rcout << "line 499. \n" ;
  
  
  NumericMatrix losstemp = Losses(0);
  // Rcout << "line 503. \n" ;
  
  int Trows = losstemp.nrow();
  int Hcols = losstemp.ncol();
  // Rcout << "line 507. \n" ;
  
  int Mmethods = Losses.size();
  int Mmethods2 = Losses.size();
  
  // Rcout << "line 510. \n" ;
  
  arma::cube Set(Trows, Hcols, Mmethods );
  
  for(int i=0;i<Losses.size();i++){
    // NumericMatrix tempmat0 = Losses(i);
    // 
    // NumericMatrix tempmat = Rcpp::clone(tempmat0);
    Set.slice(i) = as<arma::mat>(Losses(i));

  }
  
  // Rcout << "line 517. \n" ;
  
  arma::vec p_values = arma::zeros<arma::vec>(Mmethods);
  
  arma::uvec IDs = arma::regspace<arma::uvec>(0,Mmethods-1); 
  
  //This part of the code should be written in an R wrapper function
  // if(unif_or_avg =="a" &  is.null(weights)){
  //   weights <- rep(1/Hcols, Hcols)
  // }
  
  //naming uSPA, but can be aSPA depending on 
  arma::mat  t_uSPA = arma::zeros<arma::mat>(Mmethods,Mmethods);
  arma::mat  c_uSPA = arma::zeros<arma::mat>(Mmethods,Mmethods);
  
  //eventually parallelize this using number of cores equal to the minimum of the number available and the number of methods
  
  
  // std::random_device device;
  // dqrng::xoshiro256plus gen(device()); 
  
  dqrng::xoshiro256plus gen(seed); 
  
  
  
  std::uniform_real_distribution<> dis_cont_unif(0, 1);
  
  // Rcout << "line 548. \n" ;
  
  //eventually parallelize this using number of cores equal to the minimum of the number available and the number of methods
  int ncoretemp = std::min(ncores, Mmethods);
  
  // Rcout << "ncoretemp = " << ncoretemp << ". \n" ;
  
#pragma omp parallel num_threads(ncoretemp)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps
  
  // Rcout << "line 552. \n" ;
  
  
#pragma omp for
  for(int i=0; i < Mmethods;i++){
    // Rcout << "line 557. \n" ;
    
    for(int j=0; j < Mmethods;j++){
      if(i!=j){
        arma::mat LossDiff = Set.slice(i) - Set.slice(j);
        // Rcout << "line 562. \n" ;
        
        if(unif_or_avg==1){ //1 means uniform
          
          //eventually will replace this with all the function
          //to ensure use rng properly
          // arma::field<arma::vec> bootout = Bootstrap_uSPA_cpp(LossDiff, L, B) ; 
          
          
          // Rcout << "line 294 . \n";
          
          //vector of column means
          arma::rowvec d_ij = arma::mean(LossDiff);
          
          // Rcout << "line 296 . \n";
          
          //should output a vector of length 1
          //numerator is vector of length H (num columns)
          //denominator is vector of length H (num columns)
          // element-wise division, then minimum over vector
          
          //perhaps arma::min returns a double
          //arma::vec t_uSPA = arma::min(std::sqrt(TT)*d_ij/arma::sqrt(QS_cpp(Loss_Diff))) ;
          
          //arma::vec t_uSPA(1);
          t_uSPA(i,j) = arma::min(std::sqrt(Trows)*d_ij.t()/arma::sqrt(QS_cpp(LossDiff))) ;
          
          arma::vec tempsamps_t_uSPA_b = arma::zeros<arma::vec>(B);
          
          
          // Rcout << "line 315 . \n";
          
          arma::mat Demeaned_Loss_Diff = LossDiff.each_row() - d_ij;
          
          // Rcout << "line 591. \n" ;
          
          for(int b=0; b<B;b++){
            // arma::uvec id = Get_MBB_ID_cpp(TT,L);
            
            arma::uvec id = arma::zeros<arma::uvec>(Trows);
            
            
            // std::random_device device;
            // 
            // dqrng::xoshiro256plus gen(device());              // properly seeded rng
            //dqrng::xoshiro256plus gen(seed);              // properly seeded rng
            
            
            double temprand = dis_cont_unif(lgen);
            
            id(0) = std::ceil((double)Trows*temprand)-1;
            
            for(int t=1; t<Trows;t++){
              if( (t+1) % L == 0){
                id(t) = std::ceil(Trows*dis_cont_unif(lgen))-1;
              }else{
                id(t) = id(t-1);
              }
              
              if(id(t) > Trows-1){
                id(t) = 0 ;
              }
              
            }//end of for loop over t
            
            // Rcout << "line 350 . \n";
            
            arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id);
            arma::vec omega_b = MBB_Variance_cpp(b_lossdiff, L);
            
            // Rcout << "line 355 . \n";
            
            tempsamps_t_uSPA_b(b) = arma::min(std::sqrt(Trows)*( arma::mean(b_lossdiff).t())/arma::sqrt(omega_b));
            
            
          }
          
          // Rcout << "line 359 . \n";
          
          
          // Rcout << "line 637. \n" ;
          
          // vector with one element
          // arma::vec tempvec = bootout(0);
          // t_uSPA(i,j) = tempvec(0);
          // arma::vec tempsamps = bootout(1);
          
          //outputs vector of length 1
          arma::vec tempquantvec = arma::quantile(tempsamps_t_uSPA_b , alpha_t_vec) ;
          
          c_uSPA(i,j) = tempquantvec(0);
        }else{//else average
          // arma::field<arma::vec> bootout = Bootstrap_aSPA_cpp(LossDiff, weights_arma, L, B) ; 
          
          
          //really a vector, saving as matrix for similarity to R code
          arma::mat Weighted_Loss_Diff = arma::sum((LossDiff.each_row() % weights_arma.t()),1);
          
          //instead of column means, appears to be one mean of weighted row sums
          double d_ij = arma::as_scalar((arma::mean(Weighted_Loss_Diff)));
          
          
          
          //vector of length one. Keeping as vector for field output
          arma::vec t_aSPA = std::sqrt(Trows)*d_ij/arma::sqrt(QS_cpp(Weighted_Loss_Diff));
          
          t_uSPA(i,j) = t_aSPA(0);
          
          
          
          //double t_aSPA = arma::as_scalar(std::sqrt(Trows)*d_ij/arma::sqrt(QS_cpp(Weighted_Loss_Diff)));
          
          arma::vec tempsamps_t_aSPA_b = arma::zeros<arma::vec>(B);
          
          // demean the vector (saved as matrix)
          arma::mat Demeaned_Loss_Diff = Weighted_Loss_Diff - d_ij;
          
          // Rcout << "line 687. \n" ;
          
          for(int b=0; b<B;b++){
            //arma::uvec id = Get_MBB_ID_cpp(Trows,L );
            
            arma::uvec id = arma::zeros<arma::uvec>(Trows);
            
            
            // std::random_device device;
            // 
            // dqrng::xoshiro256plus gen(device());              // properly seeded rng
            //dqrng::xoshiro256plus gen(seed);              // properly seeded rng
            
            double temprand = dis_cont_unif(lgen);
            
            id(0) = std::ceil((double)Trows*temprand)-1;
            
            for(int t=1; t<Trows;t++){
              if( (t+1) % L == 0){
                id(t) = std::ceil(Trows*dis_cont_unif(lgen))-1;
              }else{
                id(t) = id(t-1)+1;
              }
              
              if(id(t) > Trows-1){
                id(t) = 0 ;
              }
              
            }//end of for loop over t
            
            // Rcout << "id = " << id << ". \n";
            
            arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id);
            
            // Rcout << "line 260 . \n";
            
            arma::vec zeta_b = MBB_Variance_cpp(b_lossdiff, L);
            //only one column, so arma::mean gives vector of length one
            tempsamps_t_aSPA_b(b) = arma::as_scalar(std::sqrt(Trows)*( arma::mean(b_lossdiff))/arma::sqrt(zeta_b));
          }
          
          
          // // vector with one element
          // arma::vec tempvec = bootout(0);
          // t_uSPA(i,j) = tempvec(0);
          // arma::vec tempsamps = bootout(1);
          
          //outputs vector of length 1
          arma::vec tempquantvec = arma::quantile(tempsamps_t_aSPA_b , alpha_t_vec) ;
          
          c_uSPA(i,j) = tempquantvec(0);
        }//end else statement
        // Rcout << "line 725. \n" ;
        
        
      }//end if statement
    }//end loop over j
  }// end loop over i
  
}//end of pragma omp code
  
      
      // Rcout << "line 740. \n" ;

      
  arma::cube t_uSPA_b = arma::zeros<arma::cube>(Mmethods, Mmethods, B );
  arma::cube c_uSPA_b = arma::zeros<arma::cube>(Mmethods, Mmethods, B );
  
  
  //eventually parallelize this using number of cores equal to the minimum of the number available and B
  ncoretemp = std::min(ncores, B);
  
   // Rcout << "line 758. \n" ;
  
  
#pragma omp parallel num_threads(ncoretemp)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps
  // 
#pragma omp for
  for(int b=0; b < B;b++){
     // Rcout << "line 768. \n" ;
    
    arma::uvec id = arma::zeros<arma::uvec>(Trows);
    
    
    // std::random_device device;
    // 
    // dqrng::xoshiro256plus gen(device());              // properly seeded rng
    //dqrng::xoshiro256plus gen(seed);              // properly seeded rng
    
    
    double temprand = dis_cont_unif(lgen);
    
    id(0) = std::ceil((double)Trows*temprand)-1;
    
    for(int t=1; t<Trows;t++){
      if( (t+1) % L == 0){
        id(t) = std::ceil(Trows*dis_cont_unif(lgen))-1;
      }else{
        id(t) = id(t-1);
      }
      
      if(id(t) > Trows-1){
        id(t) = 0 ;
      }
      
      
    }//end of for loop over t
    
    // Rcout << "line 797. \n" ;
    
    ///////////////////////////////////////////////
    //obtained IDs
    /////////////////////////
    
    for(int i=0; i < Mmethods;i++){
      for(int j=0; j < Mmethods;j++){
        if(i!=j){
          
          arma::mat LossDiff = Set.slice(i) - Set.slice(j);
          // Rcout << "line 808. \n" ;
          
          int ncoltemp = 1;
          
          if(unif_or_avg == 0){ //average
            ncoltemp = Hcols;
          }
          
          arma::mat Demeaned_Loss_Diff(Trows, ncoltemp );
          
          double t_uSPA_temp = 1;
          if(unif_or_avg == 0){ //average
            
            // arma::mat LossDiff = Set.slice(i) - Set.slice(j);
            LossDiff = arma::sum((LossDiff.each_row() % weights_arma.t()),1);
            
            // arma::mat Demeaned_Loss_Diff = LossDiff - arma::as_scalar(mean(LossDiff));
            
            Demeaned_Loss_Diff = (LossDiff - arma::as_scalar(arma::mean(LossDiff)));
            Demeaned_Loss_Diff = Demeaned_Loss_Diff.rows(id);

            // double d_ij = arma::as_scalar(mean(Demeaned_Loss_Diff));
            double d_ij = arma::as_scalar((arma::mean(Demeaned_Loss_Diff)));
            
            // t_uSPA_temp = arma::as_scalar( std::sqrt(Trows)*d_ij/arma::sqrt(QS_cpp(Demeaned_Loss_Diff)) );
            t_uSPA_temp = arma::as_scalar( std::sqrt(Trows)*d_ij/arma::sqrt(QS_cpp(Demeaned_Loss_Diff)) );
            
            
          }else{ //else uniform
            // arma::mat LossDiff = Set.slice(i) - Set.slice(j);
            
            // Rcout << "line 821. \n" ;
            
            
            // arma::mat Demeaned_Loss_Diff = LossDiff.each_row() - arma::mean(LossDiff);
            Demeaned_Loss_Diff = (LossDiff.each_row() - arma::mean(LossDiff));
            Demeaned_Loss_Diff = Demeaned_Loss_Diff.rows(id);
            
            // Rcout << "line 828. \n" ;
            // Rcout << "Demeaned_Loss_Diff.n_rows " << Demeaned_Loss_Diff.n_rows << " \n" ;
            
            
            // arma::rowvec d_ij = arma::mean(Demeaned_Loss_Diff);
            arma::rowvec d_ij = arma::mean(Demeaned_Loss_Diff);
            
            //arma::vec t_uSPA(1);
            //t_uSPA(0) =
            // Rcout << "line 850. \n" ;
            
            // t_uSPA_temp  = arma::min(std::sqrt(Trows)*d_ij.t()/arma::sqrt(QS_cpp(Demeaned_Loss_Diff))) ;
            t_uSPA_temp  = arma::min(std::sqrt(Trows)*d_ij.t()/arma::sqrt(QS_cpp(Demeaned_Loss_Diff))) ;
            
            
          }//end else statement
          
          // Rcout << "line 851. \n" ;
          
          //now apply bootstrap again
          
          arma::vec t_uSPA_b_temp = arma::zeros<arma::vec>(B);
          
          for(int b2=0; b2 < B;b2++){
            arma::uvec id2 = arma::zeros<arma::uvec>(Trows);
            
            double temprand = dis_cont_unif(lgen);
            
            id2(0) = std::ceil((double)Trows*temprand)-1;
            
            for(int t=1; t<Trows;t++){
              if( (t+1) % L == 0){
                id2(t) = std::ceil(Trows*dis_cont_unif(lgen))-1;
              }else{
                id2(t) = id2(t-1);
              }
              
              if(id2(t) > Trows-1){
                id2(t) = 0 ;
              }
            }//end of for loop over t
            
            // arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id2);
            arma::mat b_lossdiff = Demeaned_Loss_Diff.rows(id2);
            
            // Rcout << "line 879. \n" ;
            
            
            arma::vec omega_b = MBB_Variance_cpp(b_lossdiff, L);
            
            if(unif_or_avg == 0){ //average
              t_uSPA_b_temp(b2) = arma::as_scalar(std::sqrt(Trows)*( arma::mean(b_lossdiff))/arma::sqrt(omega_b));
            }else{  //uniform
              t_uSPA_b_temp(b2) = arma::min(std::sqrt(Trows)*( arma::mean(b_lossdiff).t())/arma::sqrt(omega_b));
            }
            
          }//end of inner bootsyrap over b2
          
          // Rcout << "line 892. \n" ;
          
          
          //save to relevant element of cube in outer bootstrap
          t_uSPA_b(i,j,b) = t_uSPA_temp;
          
          arma::vec tempquantvec = arma::quantile(t_uSPA_b_temp , alpha_t_vec) ;
          
          c_uSPA_b(i,j,b) = tempquantvec(0);
  
        }//end of if statement i != j
      }//end of loop over j
    }//end of loop over i
  }//end loop over b
  
}//end of pragma omp code

 // Rcout << "line 909. \n" ;
  
  
  
  while(Mmethods > 0){
    double t_max_uSPA = arma::as_scalar(arma::max(arma::max(t_uSPA	- c_uSPA)));
    
    // arma::mat t_max_uSPA_b = arma::max(t_uSPA_b - c_uSPA_b, 2);
    
    // arma::mat testmat = arma::max(x, 0);
    
    //maximum over rows, then over new rows (which were columns)
    //this gives maximum of the matrix obtained for each bootstrap sample
    arma::vec t_max_uSPA_b = arma::max(arma::max(t_uSPA_b - c_uSPA_b, 0), 0);
    
    
    arma::mat crits = t_uSPA	- c_uSPA - 9999999*arma::eye(Mmethods,Mmethods);
    
    arma::uword index =  arma::index_max(arma::max( crits , 1) );
    
    //might need to debug above line if multiple maximums (or NA values?)
    //if(length(index)==0){index <- 1} #replicating code. Perhaps should throw error here?

    arma::uvec temp_ineqs = ( t_max_uSPA < t_max_uSPA_b ) ;
    arma::vec  temp_ineqs2 = arma::conv_to<arma::vec>::from(temp_ineqs);
    
    double p_temp = arma::mean( temp_ineqs2 );
    
    double temp_max = arma::max(p_values); 
    // arma::join_cols(p_temp.t() , p_values) ;
    
    p_values(IDs(index))  = std::max(p_temp, temp_max);
    
    //p_values(IDs(index))  = std::max(temp_pvals);
    
    
    if(IDs.n_elem==1){
      p_values(IDs(index))  = 1;
    }
    
    IDs.shed_row(index);
    
    t_uSPA.shed_row(index);
    t_uSPA.shed_col(index);
    
    c_uSPA.shed_row(index);
    c_uSPA.shed_col(index);
    
    t_uSPA_b.shed_row(index);
    t_uSPA_b.shed_col(index);
    
    c_uSPA_b.shed_row(index);
    c_uSPA_b.shed_col(index);
    
    Mmethods = Mmethods -1 ;
    
    
  } // end of while loop
  
  
  
  arma::uvec IDs2 = arma::regspace<arma::uvec>(1,Mmethods2); 
  arma::vec  temp_ids = arma::conv_to<arma::vec>::from(IDs2);
  
  
  arma::mat MCS_set_temp = arma::join_rows( temp_ids , p_values) ;
  
  arma::mat MCS_set = MCS_set_temp.rows( arma::find(p_values >= 1-alpha_mcs));
        
  //MCS_set <- MCS_set[ !( p_values < 1-alpha_mcs)  ,]
        
  List ret(2);
  
  ret[0] = wrap(MCS_set);
  ret[1] = wrap(p_values);
  
  
  return(ret);
  
}




// 
// //######################################################################################################################//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::vec get_original_arma(double low,double high,double sp_low,double sp_high,arma::vec sum_preds){
//   arma::vec original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
//   
//   return(original_y);
// }
// 
// 
// //######################################################################################################################//
// 
// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// // [[Rcpp::export]]
// arma::mat arma_sub(arma::mat x, arma::uvec pos) {
//   
//   arma::mat submat = x.rows(pos) ;  // Subset by element position and set equal to
//   
//   return submat;
// }
// 
