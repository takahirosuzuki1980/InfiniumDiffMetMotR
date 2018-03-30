MotScreeningJC <- function(infile="sel_processed_Mval.txt", outname="screening_result", ControlColnum, TreatmentColnum, MethylDemethyl="Demethyl"){
	cat("Motif list construction......\n")
	library("MotifDb")
	targetDB <- "JASPAR_CORE"
	targetORG <- c("Hsapiens", "Mmusculus")
	motifDB <- query(MotifDb, targetDB)        #extraction of motif list of "JASPER_CORE"
	motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))        #extraction of motifs of "Hsapiens" and "Mmusclus"		
	motifDBList <- as.list(motifDB)
	motGene <- values(motifDB)[,4]

	cat("Reading data......\n")
	selDataMatrix <- read.table (infile)
	cat("DMP identification......\n")  
	#extraction demethylated probes  
	if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
	 diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) >=2)
	}else if ((MethylDemethyl == "Methyl" )|| (MethylDemethyl == "M")){
	 diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) <=-2) 
	}

	DMP_IDs <- rownames(selDataMatrix )[diff_table]
	nDMP_IDs <- length(DMP_IDs)
	allProbe_IDs <- rownames(selDataMatrix)

	target_position <- probeID2position(DMP_IDs)        #conversion of DMP IDs to position
	randomProbe_IDs <- stratSampling(DMP_IDs, allProbe_IDs)        #stratified sampling for negative control
	random_position <- probeID2position(randomProbe_IDs)        #conversion of NC probe IDs to position
	positionsList <- list("target" = target_position, "random" = random_position)    #integrate DMP positoins and NC positions

	##Read human hg19 genomic sequence
	library("BSgenome.Hsapiens.UCSC.hg19")
	tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep=":"))
	genome <- eval(parse(text=tmp))

	## sequence extraction
	cat("Retreave the sequences......\n")

	seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
	sequences <- lapply(positionsList , function(x){seqExtract(positions = x, genome = genome, seq_range)})

	##make a templrary outpput directory
	tempDir <- paste(outname, "_temp", sep="")
	dir.create(tempDir)

	##writing the sequences to splitted files
	seqs <- sequences$target
	target_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, tempDir=tempDir, output_file="target" )
	seqs <- sequences$random
	random_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, tempDir=tempDir, output_file="random" )
	rm(sequences)
	rm(seqs)
	invisible(replicate(3, gc()))

	cat("motif search...,,,\n")
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

	 distPlotFile <- paste(Sys.Date(),'_',outname,'_mot_plot.pdf', sep="") ##output file name setting
	pdf(distPlotFile)

	parameter_matrix <- NULL
	##Computation by each motif in motifDB
	for (j in 1:length(target_positionsList)){
	  indicator_2 <- paste("               ", "MOTIF   ",motGene[j],":",j, "/", length(motGene),"      ", date(), "\n", sep="") #monitering of progress (standerd output of vidualizing motif name)
	  cat(indicator_2)

	  motif_name <- motGene[j]
	  target_mot_posi <- unlist(target_positionsList[[j]])
	  ctrl_mot_posi <- unlist(random_positionsList[[j]])

	  ## histogram plot
	  cat("Step 1 / Plot Histgram (Figure 1)......\n")
	  motHist(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, motif_name=motif_name)
	  ## motif enrichment score plot

	  cat("Step 2 / Plot Enrichment Score (Figure 2)......\n")
	  enrichment_scores <- enrichScoreDist(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, motif_name=motif_name, nDMP_IDs=nDMP_IDs, plot_draw=TRUE)

	  cat("Step 3 / Creating a count.pvalue Matrix ......\n")
	  windowSize <- 100
	  motif_counts_matrix <- countMatrix(target_mot_posi = target_mot_posi, ran_motifPosi = ctrl_mot_posi, seq_range=seq_range, windowSize = windowSize, slide = 50)

	  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
	     cat("Step 4 / Visualization of FC (Figure 3)......\n")
	     FCChangePlot(motif_name = motif_name, motif_counts_matrix = motif_counts_matrix)

	     cat("Step 5 / Poisson disribution model exact test (Figure 4)......\n")
	     PoisTestPlot(motif_name = motif_name, seq_range = seq_range)
	  }
	 
	  cat("Step 6 / Extraction of significantly enriched ranges......\n")
	  significant_ranges <- SigRange(motif_counts_matrix = motif_counts_matrix, cutoff = 0.00001, windowSize = windowSize)
	  cat("Step 7 / Significance test......\n")

	  parameters <- enrichTest(significant_ranges = significant_ranges, motif_counts_matrix = motif_counts_matrix, seq_range=seq_range, enrichment_scores=enrichment_scores)
	  parameter_matrix <- rbind(parameter_matrix, parameters)
	}
	dev.off()

	motID <- values(motifDB)[,2]
	motSource <- values(motifDB)[,3]
	motOrg <- values(motifDB)[,9]
	finOut <- cbind(motGene,
	motID,
	motSource,
	motOrg,
	ntarget_hits,
	nrandom_hits,
	parameter_matrix)
	##out put file name setting
	Mad3ResultOut <- paste(Sys.Date(),'_',outname,'_mot_analysis_result.txt', sep="")
	write.table (finOut, file=Mad3ResultOut, sep="\t", quote=F, row.names=F)
}

