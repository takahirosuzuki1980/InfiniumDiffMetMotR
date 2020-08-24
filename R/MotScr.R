#' PWM overrepresentation analysis pipline for illumina methylation array
#' 
#' pipline to analyze enrichment of given motif PWMs for differentially methylated region based on Illumina infinium methylation array data.
#' 
#' @param infile input file: normalized M-value matrix
#' @param motifDBList PWM data list
#' @param cutoff cutoff for delta-M-value
#' @param p.cutoff cutoff for Welch's t-test
#' @param outname output files ID
#' @param ControlColnum column number for Control sample(s)
#' @param TreatmentColnum  column number for target sample(s)
#' @param MethylDemethyl direction of Ddifferentially methylated regions (methylated or demethylated)
#' @param version version of Infinium Methjylation array (450 or 850(orEPIC))
#' @param sampling If sampling number is indicated and DMP is more than sampling number, analysis is ran with randomly selected DMP. If FALSE, all DMPs are analyzed.
#' @param min.score The minimum score for counting a match. Can be given as a character string containing a percentage (e.g. "85%") of the highest possible score or as a single number. 
#' 
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils read.table write.table
#' @importFrom stats runif na.omit
#' @importFrom grDevices pdf dev.off
#' 
#' @return analysis summary files and plots, Rdata
#' @keywords PWM, overrepresentation analysis
#' @export

MotScr <- function(infile,
	motifDBList, cutoff = 2,
    p.cutoff = 0.001,
    outname = "screening_result",
    ControlColnum, TreatmentColnum,
    MethylDemethyl = "Demethyl",
    version = "850",
    sampling = FALSE, min.score = "90%"
    ){
    if(length(motifDBList) == 0) stop("motifDBList not found")
    if(missing(infile)) stop("infile not found. data.frame or text file name.")

    if(is.character(infile) && grepl("txt$", infile)){
        #read the M-value text file
        cat("Reading M-value data...\n")
        selDataMatrix <- read.table(infile, sep = "\t")
    }

    cat("DMP identification...\n")
    DMP_IDs <- DmpId(selDataMatrix = selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, p.cutoff = p.cutoff, cutoff = cutoff, MethylDemethyl = MethylDemethyl)
    if(length(DMP_IDs) == 0) stop("No DMP!!")

    if(!sampling == FALSE & length(DMP_IDs) >= sampling){
        ##fast option
        cat(paste0("Analysis will be run with ", sampling, "randomly selected DMPs."))
        cat("\n")
        DMP_IDs <- DMP_IDs[floor(runif(sampling, 1, length(DMP_IDs)))]
    }

    nDMP_IDs <- length(DMP_IDs)
    allProbe_IDs <- rownames(selDataMatrix)
    if(version == "450"){
        probe_annotation <- InfiniumDiffMetMotR::Methyl450anno
    } else if((version == "EPIC") || (version == "850")){
        probe_annotation <- InfiniumDiffMetMotR::EPICanno
    }

    target_position <- na.omit(probeID2position(probe_IDs = DMP_IDs, anno_info = probe_annotation)) #conversion of DMP IDs to position
    randomProbe_IDs <- stratSampling(target_IDs = DMP_IDs, anno_info = probe_annotation) #stratified sampling for negative control
    random_position <- na.omit(probeID2position(probe_IDs = randomProbe_IDs, anno_info = probe_annotation)) #conversion of NC probe IDs to position
    positionsList <- list("target" = target_position, "random" = random_position) #integrate DMP positoins and NC positions

    ##write DMP position
    DMPtOutfile <- paste0(outname, '_DMP_position.txt')
    out_target_position <- cbind(probeID = rownames(target_position), target_position)

    write.table(out_target_position, file = DMPtOutfile, quote = FALSE, row.names = FALSE, col.names = TRUE)

    ## sequence extraction
    ## Read human hg19 genomic sequence
    genome <- BSgenome.Hsapiens.UCSC.hg19
    cat("\nGetting the sequences...\n")
    seq_range <- c(-5000, 5000) #range from the CpG position to be extracted
    sequences <- lapply(positionsList, function(x){ seqExtract(positions = x, genome = genome, seq_range) })

    ##make a templrary outpput directory
    tempDir <- paste0(outname, "_temp")
    dir.create(tempDir)
    ##writing the sequences to splitted files
    seqs <- sequences$target
    target_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "target")
    seqs <- sequences$random
    random_all_filenames <- writeSplitSeq(seqs = seqs, split_num = 2500, tempDir = tempDir, output_file = "random")
    rm(sequences)
    rm(seqs)
    invisible(replicate(3, gc()))

    cat("\n")
    cat(paste0("motif serch: Total ", length(motifDBList), " motifs"))
    cat("\n\tTarget regions...\n")
    ## ((multi-fasta file(multi-seqs) x motif) x [motif number])) x [multi-fasta file number]
    target_positionsList <- splitSeqMotDist(filenames = target_all_filenames, motif_list = motifDBList, min.score = min.score)
    ntarget_hits <- lapply(target_positionsList, function(x){ length(unlist(x)) })
    file.remove(target_all_filenames)
    cat("\tbackground regions...\n")
    random_positionsList <- splitSeqMotDist(filenames = random_all_filenames, motif_list = motifDBList, min.score = min.score)
    nrandom_hits <- lapply(random_positionsList, function(x){ length(unlist(x)) })
    file.remove(random_all_filenames)
    gc()

    cat("\nPlotting the results...\n")
    ##Plot output setting
    distPlotFile <- paste0(outname, '_plot.pdf') ##output file name setting
    pdf(distPlotFile)
    sigPlotDir <- paste0(outname, '_sig_plots')
    ifelse(!dir.exists(sigPlotDir), dir.create(sigPlotDir), FALSE) # make a directory for significantly enriched motifs
    parameter_matrix <- NULL
    All_motif_names <- names(motifDBList)

    ##Computation by each motif in motifDB
    for(j in 1:length(target_positionsList)){
        ##monitering of progress (standerd output of vidualizing motif name)
        indicator_2 <- paste("    ", "MOTIF   ", All_motif_names[j], ": ", j, "/", length(All_motif_names), "      ", date(), "\n", sep = "")
        cat(indicator_2)

        motif_name <- All_motif_names[j]
        target_mot_posi <- unlist(target_positionsList[[j]])
        ctrl_mot_posi <- unlist(random_positionsList[[j]])

        ## histogram plot
        cat("    Step 1 / Plot Histgram (Figure 1)...\n")
        motHist(target_mot_posi, ctrl_mot_posi, seq_range = seq_range, motif_name = motif_name)

        ## motif enrichment score plot
        cat("    Step 2 / Plot Enrichment Score (Figure 2)...\n")
        enrichment_scores <- enrichScoreDist(target_mot_posi, ctrl_mot_posi, seq_range = seq_range, motif_name = motif_name, nDMP_IDs = nDMP_IDs, outname = outname, plot_draw = TRUE)

        cat("    Step 3 / Creating a count.pvalue Matrix...\n")
        windowSize <- 100
        motif_counts_matrix <- countMatrix(target_mot_posi = target_mot_posi, ctrl_mot_posi = ctrl_mot_posi, seq_range = seq_range, windowSize = windowSize, slide = 50)

        if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
            cat("    Step 4 / Visualization of FC (Figure 3)...\n")
            FCChangePlot(motif_name = motif_name, motif_counts_matrix = motif_counts_matrix)

            cat("    Step 5 / Poisson disribution model exact test (Figure 4)...\n")
            PoisTestPlot(motif_counts_matrix = motif_counts_matrix, motif_name = motif_name, seq_range = seq_range)
        }

        cat("    Step 6 / Extraction of significantly enriched ranges...\n")
        significant_ranges <- SigRange(motif_counts_matrix = motif_counts_matrix, cutoff = 0.00001, windowSize = windowSize)

        cat("    Step 7 / Significance test...\n")
        parameters <- enrichTest(significant_ranges = significant_ranges, motif_counts_matrix = motif_counts_matrix, seq_range = seq_range, enrichment_scores = enrichment_scores, windowSize = windowSize)
        ifelse(!parameters["peak.Test"] == "significant" & file.exists(paste(outname, '_sig_plots/', motif_name, '.pdf', sep = "")),
        file.remove(paste0(outname, '_sig_plots/', motif_name, '.pdf')), FALSE) ##remove the sig_plot_file if enrichment is not significant
        parameter_matrix <- rbind(parameter_matrix, parameters)
    }
    dev.off()

    cat("\nResult table writing...\n")
    finOut <- cbind(names(motifDBList), ntarget_hits, nrandom_hits, parameter_matrix)
    ##out put file name setting
    ResultOut <- paste0(outname, '_mot_analysis_result.txt')
    ResultOutR <- paste0(outname, '_result.RData')
    write.table(finOut, file = ResultOut, sep = "\t", quote = F, row.names = F)
    save(nDMP_IDs, seq_range, outname, All_motif_names, positionsList, target_positionsList, random_positionsList, file = ResultOutR)
    file.remove(tempDir)
    cat("Completed!!\n")
}
