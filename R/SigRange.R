#' Extraction of significantly enriched ranges
#' 
#' Extraction of significantly enriched ranges based on poisson test p-value. The positions where PWM is significantly enriched are marged as a region.
#' 
#' @param motif_counts_matrix motif counts matrix output from countMatrix
#' @param cutoff cutoff p-value
#' @param windowSize windowsize used for p-value coumutation (should be same as windowSize parameter of countMatrix)
#' 
#' @return range of the regions where PWM is significantly enrihed.
#' 
#' @keywords Poisson test, p-value, region
#' @export

SigRange <- function(motif_counts_matrix, cutoff = 0.00001, windowSize){
   sigClass <- names(which(motif_counts_matrix[,4] < as.numeric(cutoff)))    ##significant Class (set p-value)
   posiClass <- names(which((motif_counts_matrix[,1] - motif_counts_matrix[,2]) >= 0))    ##to extract only "enriched", select only the class of target is more than ctrl
   sigPosiClass <- sigClass[sigClass %in% posiClass]    ##significant and posi Classes
   signum <- 1
   sigrangeStr <- NULL
   sigrangeEnd <- NULL
   significant_ranges <- NULL
   while(length(sigPosiClass) >= signum){
      sigrangeStr <- c(sigrangeStr, as.numeric(sigPosiClass[signum])-(windowSize/2))
      if(!is.na(sigPosiClass[(signum+1)])){    #if the next class is also "significant", combine
         while(abs((as.numeric(sigPosiClass[signum]) - as.numeric(sigPosiClass[(signum + 1)]))) <= windowSize){    #If next significant class is less than window size, combine as one region
            signum <- signum + 1
            if(signum == length(sigPosiClass)) break
         }
      }
      sigrangeEnd <- c(sigrangeEnd, as.numeric(sigPosiClass[signum]) + (windowSize/2))
      signum <- signum + 1
   }
   significant_ranges <- cbind(sigrangeStr, sigrangeEnd)    #Significantly enriched regions
   return(significant_ranges)
}
