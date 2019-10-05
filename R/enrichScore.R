#' Computation of a enrichment score
#'
#' Mixture distribution of karnel densities per region. 
#' karnel = gaussian
#' 
#' @param x A position from the DMP
#' @param y PWM distribution at a DMP region 
#' @param z motif distribution of s randomly selected probe region
#' @param nProbeList number of DMP
#' 
#' @return a enrichment score
#' 
#' @keywords enrichmen tscore, PWM, mixture distribution karnel density

enrichScore <- function(x, y, z, nProbeList){
  demethy_length <- length(y)
  random_length <- length(z)
  b_width <-50						 #bw.nrd0(integ)
  gaus_demethy <- (exp(-((((y-x)/b_width)^2)/2)))/(sqrt(2*pi))
  gaus_random <- (exp(-((((z-x)/b_width)^2)/2)))/(sqrt(2*pi))
  mixed_gaus <- sum(gaus_demethy)-sum(gaus_random)    #Background subtraction
  enrichment_score <- mixed_gaus/nProbeList    ##normalized by region number
  return(enrichment_score)
}
