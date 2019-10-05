#' conversion of probe ID to genomic corrdinate
#' 
#' #conversion from illumina methylation array ID to genomic corrdinate (hg19)
#' 
#' @param probe_IDs probe IDs of methylation array 
#' @param anno_info annotation data base (Methyl450anno for Methyl450 array and EPICanno for EPIC array)
#' 
#' @return genomic coordinates
#' 
#' @keywords probe ID, genoic coordinate, genomic position
#' @export

probeID2position <- function (probe_IDs, anno_info = InfiniumDiffMetMotR::EPICanno){
  ##conversion of target probe IDs to positions
  probes <- probe_IDs
  probes_anno <- anno_info[(anno_info[,1] %in% probes), ]
  CHR37 <- sapply(probes, function(x){probes_anno[probes_anno[,1] %in% x, "CHR"]})
  CHR37 <- paste("chr", CHR37, sep="")
  CPG37 <- sapply(probes, function(x){probes_anno[probes_anno[,1] %in% x, "MAPINFO"]})
  positions <- cbind(CHR37, CPG37)
  rownames(positions) <- probes
  return(positions)
}
