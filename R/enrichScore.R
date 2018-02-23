enrichScore <- function(x = point ,y=target,z=ctrl, nProbeList = nDMP_IDs){ #computation of enrichment score
  ##Mixture distribution of karnel densities per region
  ##karnel = gaussian
  ##x: x value of this function
  ##target: taget distribution data  to be tested
  ##ctrl: control distribution data(motif distribution of ramdomly selected seqs)
  demethy_length <- length(y)
  random_length <- length(z)
  b_width <-50						 #bw.nrd0(integ)
  gaus_demethy <- (exp(-((((y-x)/b_width)^2)/2)))/(sqrt(2*pi))
  gaus_random <- (exp(-((((z-x)/b_width)^2)/2)))/(sqrt(2*pi))
  mixed_gaus <- sum(gaus_demethy)-sum(gaus_random)  #Background subtraction
  enrichment_score <- mixed_gaus/nProbeList  ##normalized by region number
  return(enrichment_score)
}
