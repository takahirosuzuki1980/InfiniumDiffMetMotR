DmpId <- function(selDataMatrix=selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, p.cutoff=p.cutoff, cutoff= cutoff, MethylDemethyl=MethylDemethyl){ 
	if((length(ControlColnum) > 1)||(length(TreatmentColnum) > 1)){
		## In the case of comparison of muliple samples, run statistical test (welch t test)
		cat("\tUse Welch's T-test & dlta M\n")
		library("snow")
		cl <- makeCluster(16,type="SOCK")
		clusterExport(cl, "selDataMatrix", envir=environment())
		clusterExport(cl, "ControlColnum", envir=environment())
		clusterExport(cl, "TreatmentColnum", envir=environment())
		t.pvals <- unlist(parLapply(cl, 1:nrow(selDataMatrix), function(x){t.test(unlist(selDataMatrix[x,ControlColnum]), unlist(selDataMatrix[x,TreatmentColnum]), var.equal=F,paired=F)["p.value"]}))
		dM <- unlist(parLapply(cl, 1:nrow(selDataMatrix), function(x){mean(unlist(selDataMatrix[x,ControlColnum]))-mean(unlist(selDataMatrix[x,TreatmentColnum]))}))
		stopCluster(cl)
		t.sigIDs <- rownames(selDataMatrix)[which(t.pvals < p.cutoff)]
		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")){
			diffIDs <- rownames(selDataMatrix)[which(dM >= cutoff)]
		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diffIDs <- rownames(selDataMatrix)[which(dM <= -cutoff)]
		}
		DMP_IDs <- t.sigIDs[t.sigIDs %in% diffIDs]
	}else{
		cat("\tUse dlta M\n")
		if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {
			diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) >= cutoff)
		}else if ((MethylDemethyl == "Methyl" )||(MethylDemethyl == "M")){
			diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) <=-cutoff) 
		}
		DMP_IDs <- rownames(selDataMatrix)[diff_table]
	}
	nDMPs <- paste("\t<", length(DMP_IDs), " DMPs identified>", "\n", sep="")
	cat(nDMPs)
	return(DMP_IDs)
}
