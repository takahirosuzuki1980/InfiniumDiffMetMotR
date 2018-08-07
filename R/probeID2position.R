probeID2position <- function (probe_IDs, anno_info = Methyl450anno){ #conversion from illumina humanmethylation450 array ID to genomic position (R version 3.2 or later)
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
