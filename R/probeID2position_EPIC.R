probeID2position_EPIC <- function (probe_IDs, EPICanno = EPICanno){ #conversion from illumina humanEPIC methylation array ID to genomic position (R version 3.2 or later)
  ##conversion of target probe IDs to positions
  probes <- probe_IDs
  probes_anno <- EPICanno[(EPICanno[,1] %in% probes), ]
  CHR37 <- sapply(probes, function(x){probes_anno[probes_anno[,1] %in% x, "CHR"]})
  CHR37 <- paste("chr", CHR37, sep="")
  CPG37 <- sapply(probes, function(x){probes_anno[probes_anno[,1] %in% x, "MAPINFO"]})
  positions <- cbind(CHR37, CPG37)
  rownames(positions) <- probes
  return(positions)
}
