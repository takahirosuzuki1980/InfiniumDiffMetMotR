probeID2position <- function (probe_IDs){ #conversion from illumina human450 methylation array ID to genomic position (R version 3.2 or later)
  library("FDb.InfiniumMethylation.hg19")
  library("lumi")
  ##conversion of target probe IDs to positions
  hm450 <- get450k()
  probes <- hm450[probe_IDs]
  CHR37 <- as.vector(seqnames(probes))
  CPG37 <- as.vector(start(ranges(probes)))
  positions <- cbind(CHR37, CPG37)
  rownames(positions) <- probe_IDs
  return(positions)
}
