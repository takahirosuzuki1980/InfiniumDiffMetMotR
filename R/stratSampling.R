stratSampling <- function(target_IDs, allProbe_IDs){	#stratified sampling based on DMPs. categolies are CpG island, CpG shore and non-CGI/non-shore
  library("FDb.InfiniumMethylation.hg19")
  library("lumi")
  hm450 <- get450k()
  probes <- hm450[target_IDs]
  all_probes <- hm450[allProbe_IDs]
  ##Stratified sampling of same number of probes from ctrl probe list based on categolies of CGIs, shores, and non CGIs
  ##definition of CpG island and shore
  data(hg19.islands)
  hg19.shores <- c(flank(hg19.islands, 2000, start=TRUE),flank(hg19.islands, 2000, start=FALSE))	#definign of shores (defined as 2kb up/down stream from CGI)
  ## identification of DMPs in each categolies
  CGI.probes <- subsetByOverlaps(probes, hg19.islands)	#identification of CGI DMPs
  shore.probes <- subsetByOverlaps(probes, hg19.shores)	#identification of shore DMPs
  head(shore.probes)
  head(CGI.probes)
  shore.probes <- shore.probes[!shore.probes %in% CGI.probes] #In case of overlap to two categolies, CGI is priolity
  ## number of DMPs in each categilies
  nCGI.probes <- length(CGI.probes)	#number of CGI DMPs
  nshore.probes <- length(shore.probes)	#number of shore DMPs
  nnonCGI.probes <- length(probes)-nCGI.probes-nshore.probes	#number of non-CGI/non-shore DMPs
  ##Categolization of all probes into the categolies
  CGI.ctrlProbeList <-names(subsetByOverlaps(all_probes, hg19.islands))	#identification of all CGI probes
  shore.ctrlProbeList <- names(subsetByOverlaps(all_probes, hg19.shores)) #identification of all shore probes
  shore.ctrlProbeList <- shore.ctrlProbeList[!shore.ctrlProbeList %in%CGI.ctrlProbeList]#In case of overlap to two categolies, CGI is priolity
  nonCGI.ctrlProbeList <- setdiff(setdiff(allProbe_IDs, names(CGI.ctrlProbeList)), names(shore.ctrlProbeList))	#identification of all no-CGI/non-shore probes
  ##number of probes in each categolies
  nCGI.ctrlProbeList <- length(CGI.ctrlProbeList) #number of all CGI probes
  nshore.ctrlProbeList <- length(shore.ctrlProbeList)	#number of all shore probes
  nnonCGI.ctrlProbeList <- length (nonCGI.ctrlProbeList)	#number of all non-CGI/non-shore/probes
  ##ramdom sampling from each categolies
  CGI.runProbeList <- CGI.ctrlProbeList[floor(runif(nCGI.probes, 1 ,nCGI.ctrlProbeList))]	#random sampling from CGI probes
  shore.runProbeList <- shore.ctrlProbeList[floor(runif(nshore.probes, 1 ,nshore.ctrlProbeList))]	#random sampling from shre probes
  nonCGI.runProbeList <- nonCGI.ctrlProbeList[floor(runif(nnonCGI.probes, 1 ,nnonCGI.ctrlProbeList))] 	#random sampling from non-CGI/nin-shore probes
  runProbeList <- c(CGI.runProbeList, shore.runProbeList, nonCGI.runProbeList)
  return (runProbeList)
}
