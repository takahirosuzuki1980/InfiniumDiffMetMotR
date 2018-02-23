seqExtract<- function(positions, genome , seq_range = c(-5000,5000)){ #sequence extraction
  seq_length <- seqlengths(genome)
  positions <- positions[(positions[,1] != ""),]
  chrom <- positions[,1]
  start <- as.numeric(positions[,2])-(seq_range[2]*1)	## lower range from CpG
  end <- as.numeric(positions[,2])+(seq_range[2]*1)		## upper range from CpG
  seqRange <- cbind(chrom, start, end)
  seqRange <- seqRange[(as.numeric(seqRange[,2]) >= 0) & (as.numeric(seqRange[,3]) <= seq_length[seqRange[,1]]),]
  fasta2 <- getSeq(genome, name=seqRange[,1], start=as.numeric(seqRange[,2]), end=as.numeric(seqRange[,3]))
  names(fasta2) <- rownames(positions)
  return (fasta2)
}
