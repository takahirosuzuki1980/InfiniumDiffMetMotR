#' Draw a Histgram
#' 
#' Vidualize the PWM distribution for target and randomly selected regions by histgram.
#' 
#' @param target_mot_posi PWM positions from the DMPs
#' @param ctrl_mot_posi PWM positions from the probe of random rprobes
#' @param seq_range range to be analized from the DMPs
#' @param motif_name PWM ID
#' 
#' @importFrom graphics hist
#' 
#' @keywords histgram, PWM distribution
#' @export

motHist <- function(target_mot_posi, ctrl_mot_posi, seq_range, motif_name="motif"){
  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
    ##plot a histgram
    breaks <- 400
    target_hist_para <- hist(target_mot_posi, breaks=breaks, plot=FALSE)    #parameters for histgram in target sequences
    ctrl_hist_para <- hist(ctrl_mot_posi, breaks=breaks, plot=FALSE)    #parameters for histgram in random sequences
    hist_freqs <- c(target_hist_para$counts, ctrl_hist_para$counts)
    y.limit.max <- (max(hist_freqs) *1.1)    #y max of plot area
    main<-paste(motif_name, "(histgram)", sep="")    #main title
    hist(target_mot_posi, xlim=c(seq_range[1],seq_range[2]), ylim=c(0, y.limit.max), breaks=breaks, freq=TRUE, main=main, col = "#ff00ff40", border = NA, xlab="Distance from CpG", ylab="Frequency", cex.axis=0.7)
    hist(ctrl_mot_posi, xlim=c(seq_range[1],seq_range[2]),  ylim=c(0, y.limit.max), breaks=breaks, freq=TRUE, main="",col = "#0000ff40", border = NA,, xlab="", ylab="", add=TRUE, cex.axis=0.7)
  }
}
