#' Identificaion of differentially methylated probes
#' 
#' Compare two data (or data sets) conpute differentially methylated probes.
#' For single data comparison, delta-M value is used.
#' For data set comparison, both Welch's t-test and delta-M value are used.
#' 
#' @param selDataMatrix normalized M-value matrix
#' @param ControlColnum column number for Control sample(s)
#' @param TreatmentColnum  column number for target sample(s)
#' @param p.cutoff cut-off for Welch's t-test
#' @param cutoff cut-off for delta-M value
#' @param MethylDemethyl direction of Ddifferentially methylated regions (methylated or demethylated)
#' 
#' @importFrom snow makeCluster clusterExport parLapply stopCluster
#' @importFrom stats t.test
#' 
#' @return IDs of differentially methylated probes
#' 
#' @keywords differentially methylated probe
#' @export

DmpId <- function(selDataMatrix, ControlColnum, TreatmentColnum, p.cutoff = 0.05, cutoff = 2, MethylDemethyl = "Demethyl"){ 
	if((length(ControlColnum) > 1) && (length(TreatmentColnum) > 1)){
		## For comparison of muliple samples, run statistical test (welch t test)
		cat("\tComparison between multiple data sets...Use Welch's T-test & dlta M\n")

		cl <- makeCluster(4,type="SOCK")
		clusterExport(cl, "selDataMatrix", envir=environment())
		clusterExport(cl, "ControlColnum", envir=environment())
		clusterExport(cl, "TreatmentColnum", envir=environment())
		t.pvals <- unlist(parLapply(cl, 1:nrow(selDataMatrix), function(x){t.test(unlist(selDataMatrix[x,ControlColnum]), unlist(selDataMatrix[x,TreatmentColnum]), var.equal=F,paired=F)["p.value"]}))
		dM <- unlist(parLapply(cl, 1:nrow(selDataMatrix), function(x){mean(unlist(selDataMatrix[x,ControlColnum]))-mean(unlist(selDataMatrix[x,TreatmentColnum]))}))
		stopCluster(cl)
		
		## select only methyl or demethyl
		t.sigIDs <- rownames(selDataMatrix)[which(t.pvals < p.cutoff)]
		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")){
			diffIDs <- rownames(selDataMatrix)[which(dM >= cutoff)]
		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diffIDs <- rownames(selDataMatrix)[which(dM <= -cutoff)]
		}

		DMP_IDs <- t.sigIDs[t.sigIDs %in% diffIDs]

	}else if ((length(ControlColnum) == 1) && (length(TreatmentColnum) == 1)){
		cat("\tComparison between single data\tUse dlta M\n")

		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
			diff_table <- which((selDataMatrix[,ControlColnum] - selDataMatrix[,TreatmentColnum]) >= cutoff)

		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diff_table <- which((selDataMatrix[,ControlColnum] - selDataMatrix[,TreatmentColnum]) <=-cutoff) 
		}

		DMP_IDs <- rownames(selDataMatrix)[diff_table]

	}else if ((length(ControlColnum) == 1) && (length(TreatmentColnum) > 1)){
		cat("\tComparison between single control data and ultiple treatment data\tUse dlta M with mean M of treatment\n")

		treat_mean_M <- rowMeans(selDataMatrix[,TreatmentColnum])

		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
			diff_table <- which((selDataMatrix[,ControlColnum] - treat_mean_M ) >= cutoff)

		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diff_table <- which((selDataMatrix[,ControlColnum] - treat_mean_M ) <=-cutoff) 
		}

		DMP_IDs <- rownames(selDataMatrix)[diff_table]

	}else if ((length(ControlColnum) > 1) && (length(TreatmentColnum) == 1)){
		cat("\tComparison between multiple control data and single treatment data\tUse dlta M with mean M of control\n")

		ctrl_mean_M <- rowMeans(selDataMatrix[,,ControlColnum])

		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
			diff_table <- which((ctrl_mean_M - selDataMatrix[,TreatmentColnum] ) >= cutoff)

		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diff_table <- which((ctrl_mean_M - selDataMatrix[,TreatmentColnum] ) <=-cutoff)
		}

		DMP_IDs <- rownames(selDataMatrix)[diff_table]

	}else stop ("Comparison is not correct\n")
	
	nDMPs <- paste0("\t<", length(DMP_IDs), " DMPs were identified>", "\n")
	cat(nDMPs)
	return(DMP_IDs)
}
