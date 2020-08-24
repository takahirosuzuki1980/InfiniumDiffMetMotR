#' Parallel of motifDistFromCpG with a PWM
#' 
#' Run motifDistFromCpG, a motif search function, for one PWM in parallel for multiple sequences with multi-CPU.
#' Default CPU number is 4
#' 
#' @param motif A PWM
#' @param fasta2 A DNAStringSet sequnce
#' @param min.score The minimum score for counting a match. Can be given as a character string containing a percentage (e.g. "85%") of the highest possible score or as a single number. 
#' 
#' @importFrom snow makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @importMethodsFrom Biostrings reverseComplement width matchPWM start
#' 
#' @return A list of identified motif positions in a region.
#' @keywords PWM, motif seaech, parallel


motifDistMultiSeqs <- function(motif, fasta2, min.score){

    pfm_jaspar_fw <- motif #retreave a pfm from motifDB
    pfm_jaspar_rv <- Biostrings::reverseComplement(pfm_jaspar_fw) #complementary pfm
    half_length <- (width(fasta2[1]) - 1) / 2

    ## Use multi-CPUs to compute distance from CpG
    cl <- makeCluster(4, type = "SOCK")
    clusterEvalQ(cl, library("Biostrings"))
    clusterExport(cl, "fasta2", envir = environment())
    clusterExport(cl, "min.score", envir = environment())
    clusterExport(cl, "half_length", envir = environment())
    clusterExport(cl, "pfm_jaspar_fw", envir = environment())
    clusterExport(cl, "pfm_jaspar_rv", envir = environment())
    motifPositions_fw <- parLapply(cl, fasta2, function(x){ Biostrings::start(Biostrings::matchPWM(pfm_jaspar_fw, subject = x, min.score = min.score)) - half_length })
    motifPositions_rv <- parLapply(cl, fasta2, function(x){ Biostrings::start(Biostrings::matchPWM(pfm_jaspar_rv, subject = x, min.score = min.score)) - half_length })
    stopCluster(cl)
    motifPositions <- lapply(1:length(motifPositions_fw), function(x){
        merged <- c(motifPositions_fw[[x]], motifPositions_rv[[x]])
        return(merged)
    })
    names(motifPositions) <- names(motifPositions_fw)
    return(motifPositions) #motifPositions is a list of identified motif positions in a region.
}
