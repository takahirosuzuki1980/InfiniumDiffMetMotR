motifDistMultiSeqs <- function(motif, fasta2=fasta2){ ##Run motifDistFromCpG function for multiple sequence in pallarel with multi-CPUs
  ##motif search
  library("snow")
  library("Biostrings")

  pfm_jaspar_fw <- motif #retreave a pfm from motifDB
  pfm_jaspar_rv <- reverseComplement(pfm_jaspar_fw) #complementary pfm
  half_length <- (width(fasta2[1])-1)/2

  ## Use multi-CPUs to compute distance from CpG
  cl <- makeCluster(4,type="SOCK")
  clusterEvalQ(cl,library("Biostrings"))
  clusterExport(cl, "motifDistFromCpG", envir=environment())
  clusterExport(cl, "fasta2",envir=environment())
  clusterExport(cl, "half_length",  envir=environment())
  clusterExport(cl, "pfm_jaspar_fw",  envir=environment())
  clusterExport(cl, "pfm_jaspar_rv",  envir=environment())
  motifPositions_fw <- parLapply(cl, fasta2, function(x){motifDistFromCpG(pfm=pfm_jaspar_fw, DNAseq=x, half_length=half_length)})
  motifPositions_rv <- parLapply(cl, fasta2, function(x){motifDistFromCpG(pfm=pfm_jaspar_rv, DNAseq=x, half_length=half_length)})
  stopCluster(cl)
  motifPositions <- lapply(1:length(motifPositions_fw), function(x){
    merged <- c(motifPositions_fw[[x]], motifPositions_rv[[x]])
    return(merged)
  })
  names(motifPositions) <- names(motifPositions_fw)
  return(motifPositions) #motifPositions is a list of identified motif positions in a region.
}
