motifDistFromCpG <- function(pfm, DNAseq, half_length=half_length){ ##search pfm and output the distance from center
  library("Biostrings")
  pwm_hits <- matchPWM(pfm,  subject=DNAseq,  min.score="90%")## muchPWM provides the distance from 5'end of the sequence.
  if (! identical(start(pwm_hits) , integer(0))) {
    distsFromCpG <- start(pwm_hits) - half_length  ##Then the ditance from 5' is converted to that of from center (CpG position)
    return(distsFromCpG)
  }
}
