enrichScoreDist <- function(target_mot_posi, ctrl_mot_posi, seq_range=seq_range, plot_draw=TRUE){	#Plot of enrichment score and return enrichment scores
  ## Input objects are motif position fils of 1)target and 2)control
  if((length(target_mot_posi != 0)) && (length(ctrl_mot_posi != 0))){
    ranks <- seq(seq_range[1], seq_range[2], length=1001)
    ranksList <- as.list(ranks)

    ## Use clusters (multiple CPUs)
    library("snow")
    cl <- makeCluster(16,type="SOCK")
    clusterExport(cl, "enrichScore", envir=environment())
    clusterExport(cl, "ranks", envir=environment())
    clusterExport(cl, "target_mot_posi",  envir=environment())
    clusterExport(cl, "ctrl_mot_posi",  envir=environment())
    clusterExport(cl, "nDMP_IDs",  envir=environment())
    enrichment_scores <- parLapply(cl, ranks, function(x){enrichScore(x=x, y = target_mot_posi, z = ctrl_mot_posi, nProbeList = nDMP_IDs)})
    stopCluster(cl)
    enrichment_scores <- unlist(enrichment_scores)

    if(plot_draw==TRUE){	#plot the enrichment score
      library(ggplot2)
      library(reshape2)
      library(RColorBrewer)

      main<-paste(motGene[j], "(mixture distribution)", sep="") #main title
      data <- as.matrix(enrichment_scores)
      rownames(data) <- ranks
      colnames(data ) <- "ES"
      data.df <- melt(data)
      names(data.df) <- c("dist", "sample", "ES")

      g <- ggplot(data.df, aes(x = dist, y = ES,  group = sample, colour = sample))
      g <- g + geom_line()
      g <- g + scale_colour_brewer(palette = "Set1")
      #g <- g + ylim(c(-0.01,0.03))
      g <- g + xlab("Distance from CpG")
      g <- g + ylab("Enrichment Score")
      #g <- g + guides(fill=guide_legend(title=NULL))
      g <- g + geom_hline(yintercept=0)
      g <- g + ggtitle(paste(motGene[j], "(Enrichment Score)", sep="")))
      g <- g + theme(axis.text.x = element_text( color ="black",size=14), axis.text.y = element_text(color = "black", size=14), axis.ticks=element_line(color="black", size=.5), axis.title=element_text(color="black", size=14))
      g <- g + theme(panel.background=element_rect(fill="white", color=NA), panel.border=element_rect(fill=NA, color="black"), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
      #g <- g + theme(legend.position=c(.85,.85))
      #g <- g + theme(legend.background=element_rect(fill="white"))
      g <- g + theme ()
      plot(g)
    }
    return(enrichment_scores) #returen the result of moxiture distribution
  }
}
