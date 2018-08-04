FCChangePlot <- function(motif_name = motif_name, motif_counts_matrix = motif_counts_matrix){
   main<-paste(motif_name, "(Fold_change)", sep="") #main title
   plot(rownames(motif_counts_matrix), motif_counts_matrix[,3], main=main, ylab="Log2FC (Demethyl vs Random)", xlab="Distance from CpG", type="l", pch=20, cex.axis=0.7)
}
