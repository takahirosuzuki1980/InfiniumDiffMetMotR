#' Enrichment score
#'
#'Plot of enrichment score and return enrichment scores
#'The enrichment score is conputed based on modified mixture distribution of 
#' 
#' @param target_mot_posi PWM positions from the DMPs
#' @param ctrl_mot_posi PWM positions from the probe of random rprobes
#' @param seq_range range to be analized from the DMPs
#' @param motif_name PWM ID
#' @param nDMP_IDs number of DMPs
#' @param outname ID for output file
#' @param plot_draw ligical. if TRUE, output a plot of enrichment score.
#' 
#' @importFrom snow  makeCluster clusterExport parLapply stopCluster
#' @importFrom ggplot2 ggplot geom_line scale_colour_brewer xlab ylab geom_hline ggtitle theme element_text element_line element_rect element_blank aes
#' @importFrom graphics plot
#' @importFrom grDevices pdf dev.off
#' @importFrom reshape2 melt
#' 
#' @return enrucment score
#' @keywords enrichment score, PWM, mixture distribution
#' @export
#' 

enrichScoreDist <- function(target_mot_posi, ctrl_mot_posi, seq_range, motif_name="motif", nDMP_IDs, outname="enrichment_score", plot_draw=TRUE){

  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
    ranks <- seq(seq_range[1], seq_range[2], length=1001)
    ranksList <- as.list(ranks)

    ## Use clusters (multiple CPUs)
    cl <- makeCluster(4,type="SOCK")
    clusterExport(cl, "enrichScore", envir=environment())
    clusterExport(cl, "ranks", envir=environment())
    clusterExport(cl, "target_mot_posi",  envir=environment())
    clusterExport(cl, "ctrl_mot_posi",  envir=environment())
    clusterExport(cl, "nDMP_IDs",  envir=environment())
    enrichment_scores <- parLapply(cl, ranks, function(x){enrichScore(x=x, y = target_mot_posi, z = ctrl_mot_posi, nProbeList = nDMP_IDs)})
    stopCluster(cl)
    enrichment_scores <- unlist(enrichment_scores)

    if(plot_draw==TRUE){    #plot the enrichment score
      main<-paste(motif_name, "(mixture distribution)", sep="")    #main title
      data <- as.matrix(enrichment_scores)
      rownames(data) <- ranks
      colnames(data ) <- "ES"
      data.df <- melt(data)
      names(data.df) <- c("dist", "sample", "ES")

      g <- ggplot(data.df, aes(x = data.df$dist, y = data.df$ES,  group = data.df$sample, colour = data.df$sample))
      g <- g + geom_line()
      g <- g + scale_colour_brewer(palette = "Set1")
      #g <- g + ylim(c(-0.01,0.03))
      g <- g + xlab("Distance from CpG")
      g <- g + ylab("Enrichment Score")
      #g <- g + guides(fill=guide_legend(title=NULL))
      g <- g + geom_hline(yintercept=0)
      g <- g + ggtitle(paste(motif_name, "(Enrichment Score)", sep=""))
      g <- g + theme(axis.text.x = element_text( color ="black",size=14), axis.text.y = element_text(color = "black", size=14), axis.ticks=element_line(color="black", size=.5), axis.title=element_text(color="black", size=14))
      g <- g + theme(panel.background=element_rect(fill="white", color=NA), panel.border=element_rect(fill=NA, color="black"), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
      #g <- g + theme(legend.position=c(.85,.85))
      #g <- g + theme(legend.background=element_rect(fill="white"))
      g <- g + theme ()
      plot(g)

      sigPlotFile <- paste(outname,'_sig_plots/',motif_name,'.pdf', sep="")    #output file name setting
	    pdf(file.path(sigPlotFile))
      plot(g)
      dev.off()
    }
    return(enrichment_scores)
  }
}
