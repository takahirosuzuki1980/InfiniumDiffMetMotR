#' PWM number comparison
#' 
#' Count PWM in a window and compute fold-change(FC) and poisson p-value.
#' 
#' @param target_mot_posi PWM positions of target regions
#' @param ctrl_mot_posi PWM positions of randomely selected regions
#' @param seq_range range from probe position to be analyze
#' @param windowSize window size
#' @param slide length of window sliding
#' 
#' @importFrom stats ppois p.adjust
#' 
#' @return a matrix containg number of PWM at demethylated regions, number of PWM at random regions, TF, and adjust p-value
#' 
#' @keywords count of PWM, fold-cange, p-value
#' @export

countMatrix <- function(target_mot_posi, ctrl_mot_posi, seq_range, windowSize = 100, slide = 50){
      ## motif counting
      windowSize <- windowSize    ##widow size
      slide <- slide    ##slide
      demethyMotifCounts <- motifClassCount(motifPosi=target_mot_posi, windowSize = windowSize, slide = slide, seq_range=seq_range)    #counting of motifs
      ranMotifCounts <- motifClassCount(motifPosi=ctrl_mot_posi, windowSize = windowSize, slide = slide, seq_range=seq_range)    #counting of motifs
      motif_counts_matrix <- NULL
      motif_counts_matrix <- cbind(demethyMotifCounts[,2]+1, ranMotifCounts[,2]+1,
            ((demethyMotifCounts[,2]+1)/(ranMotifCounts[,2]+1)))    #data Matrix / to avoid lambda = 0 and FC=infnity add 1
      rownames(motif_counts_matrix) = demethyMotifCounts[,1]

      ##poisson p-value Computation
      ##TEST: Poisson distribution model based test
      poisModP_unadjust <- NULL
      poisModP_adjusted <- NULL
      for (posi in 1:nrow(motif_counts_matrix)){
            ##exact test based on pirsson distribution model
            poisModP_unadjust <- c(poisModP_unadjust, (ppois(motif_counts_matrix[posi,1], lambda=motif_counts_matrix[posi,2], lower.tail=FALSE)))    #p-value of poison distribution model (upper sided)
            poisModP_adjusted <- p.adjust(poisModP_unadjust, method="BH")    #adjusted p-value
      }
      motif_counts_matrix <- cbind(motif_counts_matrix,  poisModP_adjusted)    #add adjusted p-value to data Matrix
      colnames(motif_counts_matrix) <- c("Demethylated", "Random", "FC", "adjusted.P")
      return(motif_counts_matrix)
}
