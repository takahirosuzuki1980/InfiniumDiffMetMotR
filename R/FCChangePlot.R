#' Plot a Fold-Change plot
#' 
#' Vidualization of FC
#' 
#' @param motif_name PWM ID
#' @param motif_counts_matrix motif counts matrix output from countMatrix
#' 
#' @importFrom graphics plot
#' 
#' @keywords fold-cange, plot
#' @export

FCChangePlot <- function(motif_name = "motif", motif_counts_matrix){
   main<-paste(motif_name, "(Fold_change)", sep="")    #main title
   plot(rownames(motif_counts_matrix), motif_counts_matrix[,3], main=main, ylab="Log2FC (Demethyl vs Random)", xlab="Distance from CpG", type="l", pch=20, cex.axis=0.7)
}
