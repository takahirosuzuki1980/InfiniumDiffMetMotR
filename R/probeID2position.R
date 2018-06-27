probeID2position <- function (probe_IDs, Methyl450anno = Methyl450anno){ #conversion from illumina humanmethylation450 array ID to genomic position (R version 3.2 or later)
  ##conversion of target probe IDs to positions
  probes <- probe_IDs
  CHR37 <- as.vector (Methyl450anno[which(Methyl450anno[,1]%in% probes), "CHR"])
  CHR37 <- paste("chr", CHR37, sep="")
  CPG37 <- as.vector (Methyl450anno[which(Methyl450anno[,1] %in% probes), "MAPINFO"])
  positions <- cbind(CHR37, CPG37)
  rownames(positions) <- probes
  return(positions)
}
