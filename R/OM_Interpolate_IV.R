#' Calculate Interpolated Implied Volatility Grid
#'
#' Calculate interpolated IV grid based on current IV
#'
#' @param mnes a vector of moneyness (K/S)
#' @param vol a vector of implied volatilities
#' @param grid a vector of length two containing m,k (500,2)
#'
#' @return list of interpolated iv and moneyness
#'
#' @examples
#' See README.md in GitHub
#'
#' @references Code based on [updated 2021-06-27] Model-Free Implied Measures from Options Data (Data and Code)" found at https://www.vilkov.net/codedata.html
#'
#' @export
OM_Interpolate_IV <- function(mnes, vol, grid){
  # set the grid in terms of moneyness that we want to use to compute the
  # MFIV/ MFIS and other integrals as needed
  m <-  grid[1]  # use points -500:500, i.e. 1001 points for integral approximation
  k <-  grid[2]  # use moneyness 1/(1+2) to 1+2 - should be enough as claimed by JT
  u <-  (1 + k)^(1 / m)

  # create strikes using s=1, i.e. get k/s = moneyness
  mi <- seq(-m, m)
  ki <- u^mi
  iv <- rep(NA,length(ki))

  k_s_max <- max(mnes)  # OTM calls
  k_s_min <- min(mnes)  # OTM puts
  iv_max <- vol[1]  # for more OTM puts i.e we have iv_max for min k/s, i.e. for OTM put option
  iv_min <- vol[length(vol)]  # for more OTM calls  i.e. we have iv_min for OTM call option

  # calculate the interpolated/ extrapolated IV for these k/s
  ks_larger_ind = ki > k_s_max  # more OTM call
  ks_smaller_ind = ki < k_s_min  # more OTM put
  ks_between_ind = (ki >= k_s_min) & (ki <= k_s_max)

  if(any(ks_larger_ind)){iv[ks_larger_ind] <- iv_min}

  if(any(ks_smaller_ind)){iv[ks_smaller_ind] <- iv_max}

  # evaluate the spline at ki[ks_between_ind]
  if(any(ks_between_ind)){
    s <- signal::pchip(x=mnes, y=vol,xi=ki[ks_between_ind])
    iv[ks_between_ind] <- s
  }
  result <- list(iv=iv, ki=ki)
  return(result)
}
