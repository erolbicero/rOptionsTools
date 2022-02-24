#' Calculate Model Free Implied Skew & Kurtosis
#'
#' Calculate model free implied skew and kurtosis based on BKM (2003)
#'
#' @param m scalar number of grid points
#' @param k scalar grid limit
#' @param u scalar; (1 + k)^(1 / m)
#' @param ic vector of exponents for OTM calls
#' @param ip vector of exponents for OTM  puts
#' @param currcalls vector of OTM call prices
#' @param currputs vector of OTM put prices
#' @param er scalar; exp(mat * zero_rate)
#' @param V scalar; moments_mfiv_bkm * mat
#'
#' @return list of model free implied skew and kurtosis
#'
#' @examples
#' See README.md in GitHub
#'
#' @references Code based on [updated 2021-06-27] Model-Free Implied Measures from Options Data (Data and Code)" found at https://www.vilkov.net/codedata.html
#'
#' @export
kurt_skew <- function(...){
  UseMethod("kurt_skew")
}

#' @export
kurt_skew.data.table <- function(options_surface, zero_rate
                                 , m = quote((nrow(options_surface)-1)/k)
                                 , k = 2 #typical
                                 , u = quote((1 + k)^(1 / m))
                                 , ic = quote(seq(-m, m)[options_surface$k>=1])
                                 , ip = quote(seq(-m, m)[options_surface$k <1])
                                 , currcalls = quote(options_surface$mean_price[options_surface$k>=1])
                                 , currputs = quote(options_surface$mean_price[options_surface$k<1])
                                 , er = quote(exp(mat * zero_rate))
                                 , V = quote(calc_moments_mfiv_bkm(otmPrice=options_surface$mean_price,moneyness=options_surface$k,er,mat) * mat)
                                 ){
  functionArgs <- as.list(match.call())[-1] #first element is empty
  functionArgs <- lapply(functionArgs,eval)
  defaultArgs <- formals()
  defaultArgs <- defaultArgs[!(names(defaultArgs) %in% names(functionArgs))]
  defaultArgs <- c(functionArgs, defaultArgs)
  rm(functionArgs)

  for(i in names(defaultArgs)){
    if(class(defaultArgs[[i]])[1]=="call"){
      defaultArgs[[i]] <- eval(eval(defaultArgs[[i]]),envir = list2env(defaultArgs))
    }
  }

  defaultArgs[["options_surface"]] <- NULL
  defaultArgs[["zero_rate"]] <- NULL

  result <- do.call(kurt_skew.default,defaultArgs)

  return(result)

}

#' @export
kurt_skew.default <- function(m, k, u, ic, ip, currcalls, currputs, er, V){

  a <- 3 * (u - 1) * log(1 + k) / m
  b1 <- sum(ic * (2 - (log(1 + k) / m) * ic) * currcalls / u^ic)
  b2 <- sum(ip * (2 - (log(1 + k) / m) * ip) * currputs / u^ip)
  W <- a * (b1 + b2)

  a <- 4 * (u - 1) * (log(1 + k) / m)^2
  b1 <- sum(ic^2 * (3 - (log(1 + k) / m) * ic) * currcalls / u^ic)
  b2 <- sum(ip^2 * (3 - (log(1 + k) / m) * ip) * currputs / u^ip)
  X <- a * (b1 + b2)

  mu <- er - 1 - er / 2 * V - er / 6 * W - er / 24 * X
  c_ <- (er * V - mu^2)

  mfis <- (er * W - 3 * mu * er * V + 2 * mu^3) / (c_^(3 / 2))
  mfik <- (er * X - 4 * mu * er * W + 6 * er * mu^2 * V - 3 * mu^4) / c_^2

  result <- list(mfis=mfis,mfik=mfik)

  return(result)
}

#' @noRd
#' @keywords internal
calc_moments_mfiv_bkm <- function(otmPrice,moneyness,er,mat){
  Q <- otmPrice
  ki <- moneyness
  otmcalls <- moneyness >= 1
  otmputs <- !otmcalls
  rm(otmPrice,moneyness)

  dKi <- rep(NA,length(Q))

  x <- (ki[-c(1,2)] - ki[1:(length(ki)-2)]) / 2
  dKi[-c(1,length(dKi))] <- x
  dKi[1] <-  ki[2] - ki[1]
  dKi[length(dKi)] <-  ki[length(ki)] - ki[(length(ki)-1)]
  dKi <- abs(dKi)

  Ksq <-  ki^2
  # inputs for semivars:

  bkm_nom <-  dKi * ((1 - log(ki)) * Q)
  bkm_ingredients <- bkm_nom / Ksq
  bkm_multiplier <-  2 * er

  # semivariance
  moments_mfivu_bkm <-  bkm_multiplier * sum(bkm_ingredients[otmcalls]) / mat
  moments_mfivd_bkm <-  bkm_multiplier * sum(bkm_ingredients[otmputs]) / mat
  # total variance
  moments_mfiv_bkm <- moments_mfivu_bkm + moments_mfivd_bkm
  result <- moments_mfiv_bkm
  return(result)
}
