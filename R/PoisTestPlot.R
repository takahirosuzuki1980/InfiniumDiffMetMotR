PoisTestPlot <- function(Pvals = motif_counts_matrix[,4],  motif_name = motif_name, seq_range = seq_range){
   pois2 <- -log10(Pvals)
   pois2 <- replace (pois2, which(pois2 == Inf), 500) #In the case of adjusted P-value=0, relace infinity to 500
   poisPylim <- c(0, max(pois2)*1.1) #y max of the plot area
   poisPxlim <- seq_range #range to be plotted
   ylab <- "Log10P" #y-axis label
   main <- paste(motif_name, "(Pois Pvalue)", sep="")
   plot(rownames(motif_counts_matrix), pois2, ylim=poisPylim, xlim=poisPxlim, ylab=ylab, xlab="Distance from CpG", main=main, type="l", pch=20, cex.axis=0.7)
   par(new=T)
   ## if the adjested p-value is significant (p <= 0.00001), the plots turn to red
   sigPois2 <- pois2[which(pois2 > 5)]
   plot(as.numeric(names(sigPois2)), sigPois2, ylim=poisPylim, xlim=poisPxlim, ylab="", xlab="", main=main, col="red", pch=21,cex.axis=0.7)
}
