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
  CHR37 <- probes_anno[, "CHR"]
  CHR37 <- paste0("chr", CHR37)
  CPG37 <- as.numeric(probes_anno[, "MAPINFO"])
  positions <- data.frame(CHR37, CPG37)
  rownames(positions) <- probes_anno[, "IlmnID"]
  return(positions)
}
