% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiHorizonMCS.R
\name{MultiHorizonMCS}
\alias{MultiHorizonMCS}
\title{Multihorizon model confidence set and p-values}
\usage{
MultiHorizonMCS(
  Losses,
  alpha_t = 0.05,
  alpha_mcs = 0.05,
  weights = NULL,
  L,
  B = 999,
  unif_or_average = "u"
)
}
\arguments{
\item{Losses}{A list containing M matrices of dimension T by H. Each matrix contains the losses for a forecasting method. Rows correspond to time periods (T rows), columns correspond to forecast horizons (H columns).}

\item{alpha_t}{Alpha level of critical values for pairwise tests. See Quadvlieg (2021).}

\item{alpha_mcs}{Alpha level of the critical value for the MCS test. Denoted by alpha tilde in Quadvlieg (2021) section 2.2}

\item{weights}{the 1 x H vector of weights for the losses at different horizons. For instance \code{weights <- matlab::ones(1,20)/20}}

\item{L}{integer, the parameter for the moving block bootstrap}

\item{B}{integer, the number of bootstrap iterations. Default 999}

\item{unif_or_average}{If equals "a", then the confidence set is calculated based on average SPA. If equals "u", then the confidence set is calculated based on uniform SPA.}
}
\value{
Returns a list of length 2. The elements of the list are:
\item{MCS_set}{A M by 2 matrix. The first column contains the indices of the methods for which the p-value is greater than or equal to 1-alpha_mcs. The second column contains the p-values that are greater than or equal to 1-alpha_mcs.}
\item{p_values}{A vector containing the p-values for all methods.}
}
\description{
Produces Multi-Horizon Model Confidence Set and p-values as described in section 2.2. of Quadvlieg (2021).
}
\examples{


Trow <- 20 
H <- 12
Mmethods <- 9
weights <- rep(1/H,H)

loss_list <- vector(mode = "list", length = Mmethods)

loss_list[[1]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[2]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[3]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[4]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[5]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[6]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[7]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[8]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[9]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[10]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)


MultiHorizonMCS(loss_list, L=3,B=5,unif_or_average = 'u')


}
