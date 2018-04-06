stratSampling_EPIC_3 <- function(target_IDs, allProbe_IDs, EPICanno=EPICanno, CpGislands=CpGislands, TSSs=TSSs){	#stratified sampling based on DMPs. categolies are CpG island, CpG shore and non-CGI/non-shore
  library("FDb.InfiniumMethylation.hg19")

  selEPICanno <- cbind(EPICanno[,c("CHR","MAPINFO","MAPINFO", "IlmnID")],
                       score=rep(0, nrow(EPICanno)),strand=rep("*", nrow(EPICanno)))
  selEPICanno[,1] <- paste("chr", selEPICanno[,1], sep="")
  colnames(selEPICanno) <-c("chr", "start", "end", "id", "score", "strand")
  selEPICanno <- na.omit(selEPICanno)
  EPICannobed <- with(selEPICanno, GRanges(chr, IRanges(start, end), strand, score, id=id))
  names(EPICannobed) <- EPICannobed$id
  probes <- EPICannobed[which(names(EPICannobed) %in% target_IDs)]
  all_probes <-EPICannobed[which(names(EPICannobed) %in% allProbe_IDs)]

  ##Stratified sampling of same number of probes from ctrl probe list based on categolies of CGIs, shores, and non CGIs
  ##definition of CpG island and shore
  CpGshores <- c(flank(CpGislands, 2000, start=TRUE),flank(CpGislands, 2000, start=FALSE))	#definign of shores (defined as 2kb up/down stream from CGI)

  ## identification of DMPs in each categolies (TSS, CpG Island, CpG Shore, others)
  TSS.probes <- names(subsetByOverlaps(probes, TSSs))  #identification of TSS DMPs
  CGI.probes <- names(subsetByOverlaps(probes[!(probes %over% TSSs)], CpGislands))	#identification of CGI DMPs within non-TSS probes
  shore.probes <- names(subsetByOverlaps(probes[!(probes %over% TSSs)], CpGshores))	#identification of shore DMPs within non-TSS probes
  shore.probes <- shore.probes[!(shore.probes %in% CGI.probes)] #Exclude overlap with CGI
  other.probes <- setdiff(setdiff(setdiff(target_IDs,  TSS.probes), CGI.probes), shore.probes)  #identification of non-TSS, non-CGI, non-shore probes

  ## number of DMPs in each categilies
  nTSS.probes <- length(TSS.probes) #number of CGI DMPs
  nCGI.probes <- length(CGI.probes)	#number of CGI DMPs
  nshore.probes <- length(shore.probes)	#number of shore DMPs
  nnonCGI.probes <- length(other.probes)	#number of non-CGI/non-shore DMPs

  ##Categolization of all probes into the categolies
  all.TSS.probes <- names(subsetByOverlaps(all_probes, TSSs))  #identification of all TSS probes
  all.CGI.probes <-names(subsetByOverlaps(all_probes[!(all_probes %over% TSSs)], CpGislands))	#identification of all CGI probes
  all.shore.probes <- names(subsetByOverlaps(all_probes[!(all_probes %over% TSSs)], CpGshores)) #identification of all shore probes
  all.shore.probes <- all.shore.probes[!(all.shore.probes %in% all.CGI.probes)]#In case of overlap to two categolies, CGI is priolity
  all.other.probes <- setdiff(setdiff(setdiff(allProbe_IDs,  all.TSS.probes), all.CGI.probes), all.shore.probes)	#identification of all non-TSS, non-CGI, non-shore probes

  ##number of probes in each categolies
  nall.TSS.probes <- length(all.TSS.probes) #number of all CGI probes
  nall.CGI.probes <- length(all.CGI.probes) #number of all CGI probes
  nall.shore.probes <- length(all.shore.probes)	#number of all shore probes
  nall.other.probes <- length (all.other.probes)	#number of all non-CGI/non-shore/probes

  ##ramdom sampling from each categolies
  run.TSS.probes <- all.TSS.probes[floor(runif(nTSS.probes, 1, nall.TSS.probes))] #TSS probes
  run.CGI.probes <- all.CGI.probes[floor(runif(nCGI.probes, 1 ,nall.CGI.probes))]	#CGI probes
  run.shore.probes <- all.shore.probes[floor(runif(nshore.probes, 1 ,nall.shore.probes))]	#shre probes
  run.other.probes <- all.other.probes[floor(runif(nnonCGI.probes, 1 ,nall.other.probes))] 	#other probes
  run.probes <- c(run.TSS.probes, run.CGI.probes, run.shore.probes, run.other.probes)
  return (run.probes)
}
