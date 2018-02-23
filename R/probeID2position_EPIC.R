probeID2position_EPIC <- function (probe_IDs){ #conversion from illumina humanEPIC methylation array ID to genomic position (R version 3.2 or later)
  ##conversion of target probe IDs to positions
  probes <- probe_IDs
  CHR37 <- as.vector (EPICanno[which(EPICanno[,1]%in% probes), "CHR"])
  CHR37 <- paste("chr", CHR37, sep="")
  CPG37 <- as.vector (EPICanno[which(EPICanno[,1] %in% probes), "MAPINFO"])
  positions <- cbind(CHR37, CPG37)
  return(positions)
}
