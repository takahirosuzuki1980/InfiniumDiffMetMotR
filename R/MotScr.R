MotScr <- function(infile="sel_processed_Mval.txt", motifDBList = motifDBList, cutoff = 2, p.cutoff = 0.001, outname="screening_result", ControlColnum, TreatmentColnum, MethylDemethyl="Demethyl", version = "450"){
	#This function is a pipline to analyze enrichment of given motif PWMs in differentially methylated probes of illumina arrays.

	cat("Reading data...\n")
	selDataMatrix <- read.table (infile)

	cat("DMP identification...\n")
	DMP_IDs <- DmpId(selDataMatrix=selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, p.cutoff=p.cutoff, cutoff= cutoff, MethylDemethyl=MethylDemethyl)
	nDMP_IDs <- length(DMP_IDs)
	allProbe_IDs <- rownames(selDataMatrix)
	if(version=="450"){
		target_position <- probeID2position(DMP_IDs)        #conversion of DMP IDs to position
		randomProbe_IDs <- stratSampling(DMP_IDs, allProbe_IDs)        #stratified sampling for negative control
		random_position <- probeID2position(randomProbe_IDs)        #conversion of NC probe IDs to position
		positionsList <- list("target" = target_position, "random" = random_position)    #integrate DMP positoins and NC positions
	}else if ((version=="EPIC")||(version=="850")){
		target_position <- probeID2position_EPIC(probe_IDs=DMP_IDs, EPICanno=EPICanno)	#conversion of DMP IDs to position
		randomProbe_IDs  <- stratSampling_EPIC(target_IDs=DMP_IDs, EPICanno=EPICanno)	#stratified sampling for negative control
		random_position<- probeID2position_EPIC(probe_IDs=randomProbe_IDs, EPICanno=EPICanno)	#conversion of NC probe IDs to position
		positionsList <- list("target" = target_position, "random" = random_position)	#integrate DMP positoins and NC positions
	}

	## sequence extraction
	##Read human hg19 genomic sequence
	library("BSgenome.Hsapiens.UCSC.hg19")
	tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep=":"))
	genome <- eval(parse(text=tmp))
	cat("Retreave the sequences...\n")
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

	cat("motif search for target regions...\n")
	## ((multi-fasta file(multi-seqs) x motif) x [motif number])) x [multi-fasta file number]
	target_positionsList <- splitSeqMotDist(filenames=target_all_filenames,  motif_list=motifDBList)
	ntarget_hits <- lapply(target_positionsList, function(x){length(unlist(x))})
	file.remove(target_all_filenames)
	cat("motif search for background regions...\n")
	random_positionsList <- splitSeqMotDist(filenames=random_all_filenames,  motif_list=motifDBList)
	nrandom_hits <- lapply(random_positionsList, function(x){length(unlist(x))})
	file.remove(random_all_filenames)
	file.remove(tempDir)
	gc()
	
	cat("Plotting the results...\n")
	##Plot output setting
	distPlotFile <- paste(outname,'_plot.pdf', sep="") ##output file name setting
	pdf(distPlotFile)
	parameter_matrix <- NULL
	All_motif_names <- names(motifDBList)

	##Computation by each motif in motifDB
	for (j in 1:length(target_positionsList)){
		##monitering of progress (standerd output of vidualizing motif name)
		indicator_2 <- paste("    ", "MOTIF   ",All_motif_names[j],": ",j, "/", length(All_motif_names),"      ", date(), "\n", sep="") 
		cat(indicator_2)

		motif_name <- All_motif_names[j]
		target_mot_posi <- unlist(target_positionsList[[j]])
		ctrl_mot_posi <- unlist(random_positionsList[[j]])

		## histogram plot
		cat("    Step 1 / Plot Histgram (Figure 1)...\n")
		motHist(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, motif_name=motif_name)
		## motif enrichment score plot

		cat("    Step 2 / Plot Enrichment Score (Figure 2)...\n")
		enrichment_scores <- enrichScoreDist(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, motif_name=motif_name, nDMP_IDs=nDMP_IDs, plot_draw=TRUE)

		cat("    Step 3 / Creating a count.pvalue Matrix...\n")
		windowSize <- 100
		motif_counts_matrix <- countMatrix(target_mot_posi = target_mot_posi, ctrl_mot_posi = ctrl_mot_posi, seq_range=seq_range, windowSize = windowSize, slide = 50)

		if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
			cat("    Step 4 / Visualization of FC (Figure 3)...\n")
			FCChangePlot(motif_name = motif_name, motif_counts_matrix = motif_counts_matrix)

			cat("    Step 5 / Poisson disribution model exact test (Figure 4)...\n")
			PoisTestPlot(motif_counts_matrix = motif_counts_matrix, motif_name = motif_name, seq_range = seq_range)
		}

		cat("    Step 6 / Extraction of significantly enriched ranges...\n")
		significant_ranges <- SigRange(motif_counts_matrix = motif_counts_matrix, cutoff = 0.00001, windowSize = windowSize)

		cat("    Step 7 / Significance test...\n")
		parameters <- enrichTest(significant_ranges = significant_ranges, motif_counts_matrix = motif_counts_matrix, seq_range=seq_range, enrichment_scores=enrichment_scores,windowSize=windowSize)
		parameter_matrix <- rbind(parameter_matrix, parameters)
	}
	dev.off()

	cat("Result table writing...\n")
	finOut <- cbind(names(motifDBList), ntarget_hits, nrandom_hits, parameter_matrix)
	##out put file name setting
	Mad3ResultOut <- paste(Sys.Date(),'_',outname,'_mot_analysis_result.txt', sep="")
	write.table (finOut, file=Mad3ResultOut, sep="\t", quote=F, row.names=F)
	save(positionsList, target_positionsList, random_positionsList, file="result.RData")
	cat("Completed!!\n")
}
