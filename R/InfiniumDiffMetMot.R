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

probeID2position_EPIC <- function (probe_IDs){ #conversion from illumina humanEPIC methylation array ID to genomic position (R version 3.2 or later)
  ##conversion of target probe IDs to positions
  probes <- probe_IDs
  CHR37 <- as.vector (EPICanno[which(EPICanno[,1]%in% probes), "CHR"])
  CHR37 <- paste("chr", CHR37, sep="")
  CPG37 <- as.vector (EPICanno[which(EPICanno[,1] %in% probes), "MAPINFO"])
  positions <- cbind(CHR37, CPG37)
  return(positions)
}

stratSampling <- function(target_IDs, allProbe_IDs){	#stratified sampling based on DMPs. categolies are CpG island, CpG shore and non-CGI/non-shore
  library("FDb.InfiniumMethylation.hg19")
  library("lumi")
  hm450 <- get450k()
  probes <- hm450[target_IDs]
  all_probes <- hm450[allProbe_IDs]
  ##Stratified sampling of same number of probes from ctrl probe list based on categolies of CGIs, shores, and non CGIs
  ##definition of CpG island and shore
  data(hg19.islands)
  hg19.shores <- c(flank(hg19.islands, 2000, start=TRUE),flank(hg19.islands, 2000, start=FALSE))	#definign of shores (defined as 2kb up/down stream from CGI)
  ## identification of DMPs in each categolies
  CGI.probes <- subsetByOverlaps(probes, hg19.islands)	#identification of CGI DMPs
  shore.probes <- subsetByOverlaps(probes, hg19.shores)	#identification of shore DMPs
  shore.probes <- shore.probes[!shore.probes %in% CGI.probes] #In case of overlap to two categolies, CGI is priolity

  ## number of DMPs in each categilies
  nCGI.probes <- length(CGI.probes)	#number of CGI DMPs
  nshore.probes <- length(shore.probes)	#number of shore DMPs
  nnonCGI.probes <- length(probes)-nCGI.probes-nshore.probes	#number of non-CGI/non-shore DMPs

  ##Categolization of all probes into the categolies
  CGI.ctrlProbeList <-names(subsetByOverlaps(all_probes, hg19.islands))	#identification of all CGI probes
  shore.ctrlProbeList <- names(subsetByOverlaps(all_probes, hg19.shores)) #identification of all shore probes
  shore.ctrlProbeList <- shore.ctrlProbeList[!shore.ctrlProbeList %in%CGI.ctrlProbeList]#In case of overlap to two categolies, CGI is priolity
  nonCGI.ctrlProbeList <- setdiff(setdiff(allProbe_IDs, names(CGI.ctrlProbeList)), names(shore.ctrlProbeList))	#identification of all no-CGI/non-shore probes

  ##number of probes in each categolies
  nCGI.ctrlProbeList <- length(CGI.ctrlProbeList) #number of all CGI probes
  nshore.ctrlProbeList <- length(shore.ctrlProbeList)	#number of all shore probes
  nnonCGI.ctrlProbeList <- length (nonCGI.ctrlProbeList)	#number of all non-CGI/non-shore/probes

  ##ramdom sampling from each categolies
  CGI.runProbeList <- CGI.ctrlProbeList[floor(runif(nCGI.probes, 1 ,nCGI.ctrlProbeList))]	#random sampling from CGI probes
  shore.runProbeList <- shore.ctrlProbeList[floor(runif(nshore.probes, 1 ,nshore.ctrlProbeList))]	#random sampling from shre probes
  nonCGI.runProbeList <- nonCGI.ctrlProbeList[floor(runif(nnonCGI.probes, 1 ,nnonCGI.ctrlProbeList))] 	#random sampling from non-CGI/nin-shore probes
  runProbeList <- c(CGI.runProbeList, shore.runProbeList, nonCGI.runProbeList)
  return (runProbeList)

}

stratSampling_EPIC <- function(target_IDs, allProbe_IDs){	#stratified sampling based on DMPs. categolies are CpG island, CpG shore and non-CGI/non-shore
  library("FDb.InfiniumMethylation.hg19")
  selEPICanno <- cbind(EPICanno[,c("CHR","MAPINFO","MAPINFO", "IlmnID")],
                       score=rep(0, nrow(EPICanno)),strand=rep("*", nrow(EPICanno)))
  selEPICanno[,1] <- paste("chr", selEPICanno[,1], sep="")
  colnames(selEPICanno) <-c("chr", "start", "end", "id", "score", "strand")
  selEPICanno <- na.omit(selEPICanno)
  EPICannobed <- with(selEPICanno, GRanges(chr, IRanges(start, end), strand, score, id=id))
  names(EPICannobed) <- EPICannobed$id

  probes <- EPICannobed[which(names(EPICannobed) %in% target_IDs)]
  all_probes <-EPICannobed[which(names(EPICannobed) %in% allProbe_IDs)]
  ##Stratified sampling of same number of probes from ctrl probe list based on categolies of CGIs, shores, and non CGIs
  ##definition of CpG island and shore
  data(hg19.islands)
  hg19.shores <- c(flank(hg19.islands, 2000, start=TRUE),flank(hg19.islands, 2000, start=FALSE))	#definign of shores (defined as 2kb up/down stream from CGI)
  ## identification of DMPs in each categolies
  CGI.probes <- subsetByOverlaps(probes, hg19.islands)	#identification of CGI DMPs
  shore.probes <- subsetByOverlaps(probes, hg19.shores)	#identification of shore DMPs
  shore.probes <- shore.probes[!shore.probes %in% CGI.probes] #In case of overlap to two categolies, CGI is priolity

  ## number of DMPs in each categilies
  nCGI.probes <- length(CGI.probes)	#number of CGI DMPs
  nshore.probes <- length(shore.probes)	#number of shore DMPs
  nnonCGI.probes <- length(probes)-nCGI.probes-nshore.probes	#number of non-CGI/non-shore DMPs

  ##Categolization of all probes into the categolies
  CGI.ctrlProbeList <-names(subsetByOverlaps(all_probes, hg19.islands))	#identification of all CGI probes
  shore.ctrlProbeList <- names(subsetByOverlaps(all_probes, hg19.shores)) #identification of all shore probes
  shore.ctrlProbeList <- shore.ctrlProbeList[!shore.ctrlProbeList %in%CGI.ctrlProbeList]#In case of overlap to two categolies, CGI is priolity
  nonCGI.ctrlProbeList <- setdiff(setdiff(allProbe_IDs, names(CGI.ctrlProbeList)), names(shore.ctrlProbeList))	#identification of all no-CGI/non-shore probes

  ##number of probes in each categolies
  nCGI.ctrlProbeList <- length(CGI.ctrlProbeList) #number of all CGI probes
  nshore.ctrlProbeList <- length(shore.ctrlProbeList)	#number of all shore probes
  nnonCGI.ctrlProbeList <- length (nonCGI.ctrlProbeList)	#number of all non-CGI/non-shore/probes

  ##ramdom sampling from each categolies
  CGI.runProbeList <- CGI.ctrlProbeList[floor(runif(nCGI.probes, 1 ,nCGI.ctrlProbeList))]	#random sampling from CGI probes
  shore.runProbeList <- shore.ctrlProbeList[floor(runif(nshore.probes, 1 ,nshore.ctrlProbeList))]	#random sampling from shre probes
  nonCGI.runProbeList <- nonCGI.ctrlProbeList[floor(runif(nnonCGI.probes, 1 ,nnonCGI.ctrlProbeList))] 	#random sampling from non-CGI/nin-shore probes
  runProbeList <- c(CGI.runProbeList, shore.runProbeList, nonCGI.runProbeList)
  return (runProbeList)
}

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

writeSplitSeq <- function(seqs, split_num, output_file){
  n_target_seq <- length(seqs)
  ref_point <- 0
  all_filenames <- NULL
  while(ref_point < n_target_seq){
    seqsplit_str <- 1 + ref_point	#start number of sequences
    seqsplit_end <-seqsplit_str + split_num -1	#end number of sequences
    if(seqsplit_end  > n_target_seq){ #for the last incomplete number set
      seqsplit_end <-n_target_seq
    }
    split_seqs <- seqs[seqsplit_str:seqsplit_end]
    filename <- paste(tempDir,"/",output_file,"_",seqsplit_str, "-", seqsplit_end,".fa.tmp" , sep="")
    all_filenames <- c(all_filenames, filename)
    writeXStringSet(split_seqs, file=filename, format="fasta")
    ref_point <- ref_point+split_num  #range sliding
  }
  return(all_filenames) #to read the written files return the file names
}

motifDistFromCpG <- function(pfm, DNAseq, half_length=half_length){ ##search pfm and output the distance from center
  library("Biostrings")
  pwm_hits <- matchPWM(pfm,  subject=DNAseq,  min.score="80%")## muchPWM provides the distance from 5'end of the sequence.
  if (! identical(start(pwm_hits) , integer(0))) {
    distsFromCpG <- start(pwm_hits) - half_length  ##Then the ditance from 5' is converted to that of from center (CpG position)
    return(distsFromCpG)
  }
}

motifDistMultiSeqs <- function(motif, fasta2=fasta2){ ##Run motifDistFromCpG function for multiple sequence in pallarel with multi-CPUs
  ##motif search
  library("snow")
  library("Biostrings")

  pfm_jaspar_fw <- motif #retreave a pfm from motifDB
  pfm_jaspar_rv <- reverseComplement(pfm_jaspar_fw) #complementary pfm
  half_length <- (width(fasta2[1])-1)/2

  ## Use multi-CPUs to compute distance from CpG
  cl <- makeCluster(16,type="SOCK")
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
  })
  names(motifPositions) <- names(motifPositions_fw)
  return(motifPositions) #motifPositions is a list of identified motif positions in a region.
}

splitSeqMotDist <- function(filenames,  motif_list){ #Splitted seqs are sequentially read and serch the pwm
  target_positionsList <- lapply(1:length(motif_list), function(i){numeric()})  #an object to save result of motif search for all motifs
  names (target_positionsList) <- names(motif_list)
  count <- 1
  for (i in filenames){ #read the multi-fasta file names one by one
    print(paste(i," is processing......", count, "/",length(filenames), sep=""))
    fasta2 <- readDNAStringSet(i)	#read the multi-fasta file as DNAStringSet
    target_posi_subList <- lapply(motif_list, function(x){motifDistMultiSeqs(motif = x, fasta2 = fasta2)})	#target_posi_subList is a list of lists of identified motif positions in a region for each motif.
    for (j in 1:length(target_posi_subList)){
      target_positionsList[[j]] <- c(target_positionsList[[j]],target_posi_subList[[j]] )	#save the result to an List object
    }
    count <- count+1
    rm(fasta2)
    print("Done...")
  }
  return(target_positionsList)
}

motHist <- function(target_mot_posi, ctrl_mot_posi, seq_range=seq_range){ #Vidualization by histgram
  ## Input objects are 1) motif position fils of target and 2)control
  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
    ##plot a histgram
    breaks <- 400
    target_hist_para <- hist(target_mot_posi, breaks=breaks, plot=FALSE) #parameters for histgram in target sequences
    ctrl_hist_para <- hist(ctrl_mot_posi, breaks=breaks, plot=FALSE)  #parameters for histgram in random sequences
    hist_freqs <- c(target_hist_para$counts, ctrl_hist_para$counts)
    y.limit.max <- (max(hist_freqs) *1.1)  #y max of plot area
    main<-paste(motGene[j], "(histgram)", sep="") #main title
    hist(target_mot_posi, xlim=c(seq_range[1],seq_range[2]), ylim=c(0, y.limit.max), breaks=breaks, freq=TRUE, main=main, col = "#ff00ff40", border = NA, xlab="Distance from CpG", ylab="Frequency", cex.axis=0.7)
    hist(ctrl_mot_posi, xlim=c(seq_range[1],seq_range[2]),  ylim=c(0, y.limit.max), breaks=breaks, freq=TRUE, main="",col = "#0000ff40", border = NA,, xlab="", ylab="", add=TRUE, cex.axis=0.7)
  }
}

enrichScore <- function(x = point ,y=target,z=ctrl, nProbeList = nDMP_IDs){ #computation of enrichment score
  ##Mixture distribution of karnel densities per region
  ##karnel = gaussian
  ##x: x value of this function
  ##target: taget distribution data  to be tested
  ##ctrl: control distribution data(motif distribution of ramdomly selected seqs)
  demethy_length <- length(y)
  random_length <- length(z)
  b_width <-50						 #bw.nrd0(integ)
  gaus_demethy <- (exp(-((((y-x)/b_width)^2)/2)))/(sqrt(2*pi))
  gaus_random <- (exp(-((((z-x)/b_width)^2)/2)))/(sqrt(2*pi))
  mixed_gaus <- sum(gaus_demethy)-sum(gaus_random)  #Background subtraction
  enrichment_score <- mixed_gaus/nProbeList  ##normalized by region number
  return(enrichment_score)

}

enrichScoreDist <- function(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, plot_draw=TRUE){	#Plot of enrichment score and return enrichment scores
  ## Input objects are motif position fils of 1)target and 2)control
  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
    ranks <- seq(seq_range[1], seq_range[2], length=1001)
    ranksList <- as.list(ranks)

    ## Use clusters (multiple CPUs)
    library("snow")
    cl <- makeCluster(16,type="SOCK")
    clusterExport(cl, "enrichScore", envir=environment())
    clusterExport(cl, "ranks", envir=environment())
    clusterExport(cl, "target_mot_posi",  envir=environment())
    clusterExport(cl, "ctrl_mot_posi",  envir=environment())
    clusterExport(cl, "nDMP_IDs",  envir=environment())
    enrichment_scores <- parLapply(cl, ranks, function(x){enrichScore(x=x, y = target_mot_posi, z = ctrl_mot_posi, nProbeList = nDMP_IDs)})
    stopCluster(cl)
    enrichment_scores <- unlist(enrichment_scores)

    if(plot_draw==TRUE){	#plot the enrichment score
      library(ggplot2)
      library(reshape2)
      library(RColorBrewer)

      yLimitMixture <- c(min(enrichment_scores*1.1), max(enrichment_scores *1.1)) #y max of the plot area
      main<-paste(motGene[j], "(mixture distribution)", sep="") #main title
      data <- as.matrix(enrichment_scores)
      rownames(data) <- ranks
      colnames(data ) <- "ES"
      data.df <- melt(data)
      names(data.df) <- c("dist", "sample", "ES")

      g <- ggplot(data.df, aes(x = dist, y = ES,  group = sample, colour = sample))
      g <- g + geom_line()
      g <- g + scale_colour_brewer(palette = "Set1")
      #g <- g + ylim(c(-0.01,0.03))
      g <- g + xlab("Distance from CpG")
      g <- g + ylab("Enrichment Score")
      #g <- g + guides(fill=guide_legend(title=NULL))
      g <- g + geom_hline(yintercept=0)
      g <- g + theme(axis.text.x = element_text( color ="black",size=14), axis.text.y = element_text(color = "black", size=14), axis.ticks=element_line(color="black", size=.5), axis.title=element_text(color="black", size=14))
      g <- g + theme(panel.background=element_rect(fill="white", color=NA), panel.border=element_rect(fill=NA, color="black"), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
      #g <- g + theme(legend.position=c(.85,.85))
      #g <- g + theme(legend.background=element_rect(fill="white"))
      g <- g + theme ()
      plot(g)
    }
    return(enrichment_scores) #returen the result of moxiture distribution
  }
}

motifClassCount <- function(motifPosi, windowSize = windowSize, slide = slide, seq_range=seq_range){
  relativePosi <- 0
  freqs <- NULL
  class <- NULL
  windowStr <- seq_range[1] + relativePosi	#intitial start position
  windowEnd <- windowStr + windowSize	#intial end position
  while(windowEnd <= seq_range[2]){	#sliding
    windowStr <- seq_range[1] + relativePosi #5' of window end
    windowEnd <- windowStr + windowSize #3' of window end
    class <- c(class, windowStr+(windowSize/2)) #class (center of window)
    freqs <- c(freqs,sum((motifPosi >= windowStr) & (motifPosi <= windowEnd))) #counting of motif frequency
    relativePosi <- relativePosi+slide #Window sliding
  }
  freqs <- freqs
  motifCounts <- cbind(class, freqs)
  return (motifCounts)
}
