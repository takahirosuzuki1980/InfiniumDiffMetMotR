motifClassCount <- function(motifPosi, windowSize = windowSize, slide = slide, seq_range=seq_range){
  relativePosi <- 0
  freqs <- NULL
  class <- NULL
  windowStr <- seq_range[1] + relativePosi	#intitial start position
  windowEnd <- windowStr + windowSize	#intial end position
  while(windowEnd <= seq_range[2]){	#sliding
    windowStr <- seq_range[1] + relativePosi #5' of window end
    windowEnd <- windowStr + windowSize #3' of window end
    class <- c(class, windowStr+(windowSize/2)) #class (center of window)
    freqs <- c(freqs,sum((motifPosi >= windowStr) & (motifPosi <= windowEnd))) #counting of motif frequency
    relativePosi <- relativePosi+slide #Window sliding
  }
  freqs <- freqs
  motifCounts <- cbind(class, freqs)
  return (motifCounts)
}
