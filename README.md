InfiniumDiffMetMotR
===================
Version: 1.0

Description: This is a R package to analyze transcription factor binding motif enrichment for differentially methylated regions.  

Last Update: 2018-3-01

Updated by : takahiro.suzuki.aa@riken.jp

Example
-------
#### 1. motif database construction  
```
library("MotifDb")
targetDB <- "JASPAR_CORE"
targetORG <- c("Hsapiens", "Mmusculus")
targetTF <- "SPI1"
motifDB <- query(MotifDb, targetDB)        #extraction of motif list of "JASPER_CORE"
motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))        #extraction of motifs of "Hsapiens" and "Mmusclus"
motifDB <- query(motifDB,targetTF)       #Extraction of motifs for target TF(s)
motifDBList <- as.list(motifDB)
```
##discription of motifs
```
motGene <- values(motifDB)[,4]
motID <- values(motifDB)[,2]
motSource <- values(motifDB)[,3]
motOrg <- values(motifDB)[,9]
```
#### 2. Identification of differentially methylated regions  
```
infile <- "sel_processed_Mval.txt"
outname <-"iPS-HPC_SPI1"
ControlColnum <- 20
TreatmentColnum <- 11
MethylDemethyl <- "Demethyl"

selDataMatrix <- read.table (infile)  
#extraction demethylated probes  
if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
 diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) >=2)
}else if ((MethylDemethyl == "Methyl" )|| (MethylDemethyl == "M")){
 diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) <=-2) 
}

DMP_IDs <- rownames(selDataMatrix )[diff_table]
nDMP_IDs <- length(DMP_IDs)
allProbe_IDs <- rownames(selDataMatrix)
```
#### 3. Extraction of DMP positions and stratified sampling  
```
target_position <- probeID2position(DMP_IDs)        #conversion of DMP IDs to position
randomProbe_IDs <- stratSampling(DMP_IDs, allProbe_IDs)        #stratified sampling for negative control
random_position <- probeID2position(randomProbe_IDs)        #conversion of NC probe IDs to position
positionsList <- list("target" = target_position, "random" = random_position)    #integrate DMP positoins and NC positions
```
### 4. Sequence extraction  
```
##Read human hg19 genomic sequence
library("BSgenome.Hsapiens.UCSC.hg19")
tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep=":"))
genome <- eval(parse(text=tmp))

## sequence extraction
cat("Retreave the sequences\n")

seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
sequences <- lapply(positionsList , function(x){seqExtract(positions = x, genome = genome, seq_range)})
```
### 5. writing the sequences  
The extracted sequences occupy much memory. To reduce the memory occupancy, the extracted sequences is written out to splitted fasta files.   
```
##make a templrary outpput directory
tempDir <- paste(outname, "_temp", sep="")
dir.create(tempDir)

##writing the sequences to splitted files
seqs <- sequences$target
target_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, output_file="target" )
seqs <- sequences$random
random_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, output_file="random" )
rm(sequences)
rm(seqs)
invisible(replicate(3, gc()))
```
### 6. motif search  
```
## ((multi-fasta file(multi-seqs) x motif) x [motif number])) x [multi-fasta file number]
target_positionsList <- splitSeqMotDist(filenames=target_all_filenames,  motif_list=motifDBList)
ntarget_hits <- lapply(target_positionsList, function(x){
	length(unlist(x))}
	)
file.remove(target_all_filenames)

random_positionsList <- splitSeqMotDist(filenames=random_all_filenames,  motif_list=motifDBList)
nrandom_hits <- lapply(random_positionsList, function(x){
	length(unlist(x))}
	)
 file.remove(random_all_filenames)
 file.remove(tempDir)
 gc()
 ```
 ### 7. vidualization & statistical analysis  
 ```
 distPlotFile <- paste(Sys.Date(),'_',outname,'_mot_plot.pdf', sep="") ##output file name setting
pdf(distPlotFile)

poispeakPosiEvals <- NULL
maxClass <- NULL
maxFC <- NULL
poisPvals <- NULL
enrichment_Q1<- NULL
enrichment_Med <- NULL
enrichment_Q3 <- NULL
enrichment_IQR <- NULL
centPeakOutL <- NULL
centPeakMax <- NULL
centPeakRatio <- NULL
peakMaxEval <- NULL

##Computation by each motif in motifDB
for (j in 1:length(target_positionsList)){
	indicator_2 <- paste("               ", "MOTIF   ",motGene[j],":",j, "/", length(motGene),"      ", date(), "\n", sep="") #monitering of progress (standerd output of vidualizing motif name)
	cat(indicator_2)

	target_mot_posi <- unlist(target_positionsList[[j]])
	ctrl_mot_posi <- unlist(random_positionsList[[j]])

	## histogram plot
	cat("Step 1 / plot histgram\n")
	motHist (target_mot_posi, ctrl_mot_posi, seq_range=seq_range)
	## motif enrichment score plot
	cat("Step 2 / plot mixture distribution\n")
	enrichment_scores <- enrichScoreDist (target_mot_posi, ctrl_mot_posi, seq_range=seq_range, plot_draw=TRUE)
	cat("Step 3 / Poisson disribution model exact test\n")
	##TEST: Poisson distribution model based test
	##motifClassCount returns frequencies of motif in a window based on windowSize and slide.

	## motif counting
		windowSize <- 100  ##widow size
	slide <- 50	##slide
	demethyMotifCounts <- motifClassCount(motifPosi=target_mot_posi, windowSize = windowSize, slide = slide, seq_range=seq_range) #counting of motifs
	ranMotifCounts <- motifClassCount(motifPosi=ctrl_mot_posi, windowSize = windowSize, slide = slide, seq_range=seq_range) #counting of motifs

	motif_counts_matrix <- NULL
	motif_counts_matrix <- cbind(demethyMotifCounts[,2]+1, ranMotifCounts[,2]+1,
	((demethyMotifCounts[,2]+1)/(ranMotifCounts[,2]+1))) #data Matrix / to avoid lambda = 0 and FC=infnity add 1
	rownames(motif_counts_matrix) = demethyMotifCounts[,1]

	##poisson p-value Computation
	poisModP_unadjust <- NULL
	poisModP_adjusted <- NULL
	for (posi in 1:nrow(motif_counts_matrix)){
		##exact test based on pirsson distribution model
		poisModP_unadjust <- c(poisModP_unadjust, (ppois(motif_counts_matrix[posi,1], lambda=motif_counts_matrix[posi,2], lower.tail=FALSE))) #p-value of poison distribution model (upper sided)
		poisModP_adjusted <- p.adjust(poisModP_unadjust, method="BH") #adjusted p-value
	}

	motif_counts_matrix <- cbind(motif_counts_matrix,  poisModP_adjusted) #add adjusted p-value to data Matrix
	colnames(motif_counts_matrix) <- c("Demethylated", "Random", "FC", "adjusted.P")

	##vidualization of the pois test
	cat("Step 4 / Visualization of the pois test\n")
	if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
	main<-paste(motGene[j], "(Fold_change)", sep="") #main title
		plot(rownames(motif_counts_matrix), motif_counts_matrix[,3], main=main, ylab="Fold-Change (Demethyl vs. random)", xlab="Distance from CpG", type="l", pch=20, cex.axis=0.7)
		pois2 <- -log10(motif_counts_matrix[,4])
		pois2 <- replace (pois2, which(pois2 == Inf), 500) #In the case of adjusted P-value=0, relace infinity to 500
		poisPylim <- c(0, max(pois2)*1.1) #y max of the plot area
		poisPxlim <- seq_range #range to be plotted
		main<-paste(motGene[j], "(Log10P)", sep="") #main title
		plot(rownames(motif_counts_matrix), pois2, ylim=poisPylim, xlim=poisPxlim, ylab=main, xlab="Distance from CpG", main="poison distribution P-value", type="l", pch=20, cex.axis=0.7)
		par(new=T)
		## if the adjested p-value is significant (p <= 0.00001), the plots turn to red
		sigPois2 <- pois2[which(pois2 > 5)]
		plot(as.numeric(names(sigPois2)), sigPois2, ylim=poisPylim, xlim=poisPxlim, ylab="", xlab="", main=main, col="red", pch=21,cex.axis=0.7)
	}

	## Extraction of significantly enriched ranges
	cat("Step 5 / Extraction of significantly enriched ranges\n")
	sigClass <- rownames(motif_counts_matrix)[which(poisModP_adjusted < 0.00001)] ##significant Class (set p-value)
	posiClass <- rownames(motif_counts_matrix)[which((motif_counts_matrix[,1] - motif_counts_matrix[,2]) >= 0)] ##to extract only "enriched", select only the class of target is more than ctrl
	sigPosiClass <- sigClass[sigClass %in% posiClass] ##significant and posi Classes
	signum <- 1
	sigrangeStr <- NULL
	sigrangeEnd <- NULL
	sigRange <- NULL
	while(length(sigPosiClass) >= signum){
		sigrangeStr <- c(sigrangeStr,as.numeric(sigPosiClass[signum])-(windowSize/2))
		if(!is.na(sigPosiClass[(signum+1)])){ # if the next class is also "significant" (to check the range)
			while(abs((as.numeric(sigPosiClass[signum]) - as.numeric(sigPosiClass[(signum+1)]))) <=100){   #100 is window size
				signum <- signum + 1
				if (signum == length(sigPosiClass)) break
			}
		}
		sigrangeEnd <- c(sigrangeEnd,as.numeric(sigPosiClass[signum])+(windowSize/2))
		signum <- signum + 1
	}
	sigRange <- cbind(sigrangeStr, sigrangeEnd)
	if(!is.null(sigRange)){		#if significant range exist
		centSigRangeInd <- which(sigRange[,1] < 0 & sigRange[,2] > 0) #extract index (row numbers) of significanrly enriched rgion which is overlapped with CpG
	}else{
		centSigRangeInd <- NULL
	}

	if(length(centSigRangeInd) != 0){
		## Max value (FC) in the center peak
		cat("Step 6 / Max value (FC) in the center peak\n")
		centSigRange <- c(sigRange[centSigRangeInd,1], sigRange[centSigRangeInd,2]) #center peak range
		rowNoMaxClass <- which(motif_counts_matrix[,3] == max(motif_counts_matrix[which(as.numeric(rownames(motif_counts_matrix)) >= centSigRange[1] & as.numeric(rownames(motif_counts_matrix)) <= centSigRange[2]),3]))  #row index of max FC
		if( length(rowNoMaxClass) > 1){ ## if the min Pvalue is not single, adopt a class which is the closest to 0
			rowNoMaxClass <- which(rownames(motif_counts_matrix) == min(rownames(motif_counts_matrix)[rowNoMaxClass]))
		}
		maxClass <- c(maxClass, rownames(motif_counts_matrix)[rowNoMaxClass])
		maxFC <- c(maxFC, motif_counts_matrix[rowNoMaxClass, 3])
		poisPvals <- c(poisPvals, min(poisModP_adjusted[which(as.numeric(rownames(motif_counts_matrix)) > sigRange[centSigRangeInd,1] & as.numeric(rownames(motif_counts_matrix)) < sigRange[centSigRangeInd,2])]))
		poispeakPosiEvals <- c(poispeakPosiEvals, "Overlapped")

		##TEST: Mixture distribution based test
		cat("Step 7 /  enrichment score center peak test\n")
		ranks <- seq(seq_range[1], seq_range[2], length=1001)
		centPeakProb <- enrichment_scores[which(ranks >centSigRange[1] & ranks < centSigRange[2])] #"enrichments" (values of mixture distribution) of center peak
		max_value <- max(centPeakProb) #max of "enrichment" in the center peak
		max_ratio <- (max_value - quantile(enrichment_scores)[4])/IQR(enrichment_scores)
		outlier <- as.numeric(quantile(enrichment_scores)[4]) + (IQR(enrichment_scores)*3) # set of outlier cut-off (Q3+IQR*3)
		enrichment_Q1<- c(enrichment_Q1, quantile(enrichment_scores)[2])
		enrichment_Med <-c(enrichment_Med, quantile(enrichment_scores)[3])
		enrichment_Q3 <- c(enrichment_Q3, quantile(enrichment_scores)[4])
		enrichment_IQR <- c(enrichment_IQR, IQR(enrichment_scores))
		centPeakRatio <- c(centPeakRatio, max_ratio)
		centPeakOutL <- c(centPeakOutL, outlier)
		centPeakMax <- c(centPeakMax, max_value)
		## main body of judgement of center peak test
		if (max_value > outlier){
			peakMaxEval <- c(peakMaxEval, "significant")
		} else {
			peakMaxEval <- c(peakMaxEval, "NS")
		}
	}else{
		##case of no center peak
		cat("Step 6 /abort: no significant center peak\n")
		poisPvals <- c(poisPvals,min(poisModP_adjusted[which(as.numeric(rownames(motif_counts_matrix)) > (0-(windowSize/2)) & as.numeric(rownames(motif_counts_matrix)) < (0+(windowSize/2)))]))  #put min P of cnter window
		poispeakPosiEvals <- c(poispeakPosiEvals, "Not-Overlapped") #Judgement of center peak
		maxClass <- c(maxClass, "NA")
		maxFC <- c(maxFC, "NA")
		centPeakMax <- c(centPeakMax, "NA")
		enrichment_Q1<- c(enrichment_Q1, "NA")
		enrichment_Med <-c(enrichment_Med, "NA")
		enrichment_Q3 <- c(enrichment_Q3, "NA")
		enrichment_IQR <- c(enrichment_IQR, "NA")
		centPeakRatio <- c(centPeakRatio, "NA")
		centPeakOutL <- c(centPeakOutL, "NA")
		peakMaxEval <- c(peakMaxEval, "NA")
	}
}
dev.off()


##Result output
finOut <- cbind(motGene,
motID,
motSource,
motOrg,
ntarget_hits,
nrandom_hits,
poispeakPosiEvals,
maxClass,
maxFC,
poisPvals,
enrichment_Q1,
enrichment_Med,
enrichment_Q3,
enrichment_IQR,
centPeakOutL,
centPeakMax,
centPeakRatio,
 peakMaxEval)
 ##out put file name setting
Mad3ResultOut <- paste(Sys.Date(),'_',outname,'_mot_analysis_result.txt', sep="")
write.table (finOut, file=Mad3ResultOut, sep="\t", quote=F, row.names=F)
```
