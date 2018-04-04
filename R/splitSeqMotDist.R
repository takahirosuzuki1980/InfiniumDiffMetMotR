splitSeqMotDist <- function(filenames,  motif_list){ #Splitted seqs are sequentially read and serch the pwm
  library("pbapply")
  target_positionsList <- lapply(1:length(motif_list), function(i){numeric()})  #an object to save result of motif search for all motifs
  names (target_positionsList) <- names(motif_list)
  count <- 1
  for (i in filenames){ #read the multi-fasta file names one by one
    sto <- paste(i," is processing......", count, "/",length(filenames), "\n", sep=""))
    cat(sto)
    fasta2 <- readDNAStringSet(i)	#read the multi-fasta file as DNAStringSet
    target_posi_subList <- pblapply(motif_list, function(x){motifDistMultiSeqs(motif = x, fasta2 = fasta2)})	#target_posi_subList is a list of lists of identified motif positions in a region for each motif.
    for (j in 1:length(target_posi_subList)){
      target_positionsList[[j]] <- c(target_positionsList[[j]],target_posi_subList[[j]] )	#save the result to an List object
    }
    count <- count+1
    rm(fasta2)
    cat("Done...\n")
  }
  return(target_positionsList)
}
