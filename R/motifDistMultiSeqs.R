#' Parallel of motifDistFromCpG with a PWM
#' 
#' Run motifDistFromCpG, a motif search function, for one PWM in parallel for multiple sequences with multi-CPU.
#' Default CPU number is 4
#' 
#' @param motif A PWM
#' @param fasta2 DNAStringSet sequnces
#' 
#' @importFrom snow makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @importFrom Biostrings reverseComplement width
#' 
#' @return A list of identified motif positions in a region.
#' 
#' @keywords PWM, motif seaech, parallel
#' @export

motifDistMultiSeqs <- function(motif, fasta2){
  ## motif search
  pfm_jaspar_fw <- motif    #retreave a pfm from motifDB
  pfm_jaspar_rv <- reverseComplement(pfm_jaspar_fw)    #complementary pfm
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
    }
  )
  names(motifPositions) <- names(motifPositions_fw)
  return(motifPositions)
}
