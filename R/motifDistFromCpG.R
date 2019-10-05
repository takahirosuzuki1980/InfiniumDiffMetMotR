#' PWM serch with singel motif and single sequence
#'
#' search a pfm and output the distance from center
#' The default minimam score = 90%
#' 
#' @param pfm A PWM
#' @param DNAseq A sequence
#' @param half_length half of sequence length
#' 
#' @importFrom Biostrings matchPWM
#' @importFrom stats start
#' 
#' @return mached PWM positions (distances from probe)
#' @keywords PWM, motif seaech
#' @export

motifDistFromCpG <- function(pfm, DNAseq, half_length){ 
  pwm_hits <- matchPWM(pfm, subject=DNAseq,  min.score="90%")    #muchPWM provides the distance from 5'end of the sequence.
  if (! identical(start(pwm_hits) , integer(0))) {
    distsFromCpG <- start(pwm_hits) - half_length    #Then the ditance from 5' is converted to that of from center (CpG position)
    return(distsFromCpG)
  }
}
