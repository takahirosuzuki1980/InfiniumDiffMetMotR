




stratSampling_EPIC <- function(target_IDs, EPICanno=EPICanno){	#stratified sampling based on DMPs. categolies are CpG island, CpG shore and non-CGI/non-shore
 ##Categolization of all probes into the categolies
  all.TSS.probes <- as.vector(EPICanno[unique(c(grep("TSS200", EPICanno[,"UCSC_RefGene_Group"]), grep("5'UTR", EPICanno[,"UCSC_RefGene_Group"]))),"IlmnID"])#identification of all TSS probes
  all.nonTSS <- EPICanno[!(EPICanno[,"IlmnID"] %in% all.TSS.probes),]
  all.CGI.probes <- as.vector(all.nonTSS[grep("Island", all.nonTSS[,"Relation_to_UCSC_CpG_Island"]),"IlmnID"])#identification of all CGI probes
  all.shore.probes <- as.vector(all.nonTSS[grep("Shore", all.nonTSS[,"Relation_to_UCSC_CpG_Island"]),"IlmnID"]) #identification of all shore probes
  all.shelf.probes <- as.vector(all.nonTSS[grep("Shelf", all.nonTSS[,"Relation_to_UCSC_CpG_Island"]),"IlmnID"]) #identification of all shelf probes
  all.other.probes <- setdiff(setdiff(setdiff(setdiff(EPICanno[,"IlmnID"],all.TSS.probes), all.CGI.probes), all.shore.probes),all.shelf.probes)  #identification of all non-TSS, non-CGI, non-shore probes

  ##number of probes in each categolies
  nall.TSS.probes <- length(all.TSS.probes) #number of all CGI probes
  nall.CGI.probes <- length(all.CGI.probes) #number of all CGI probes
  nall.shore.probes <- length(all.shore.probes) #number of all shore probes
  nall.shelf.probes <- length(all.shelf.probes) #number of all shelf probes
  nall.other.probes <- length (all.other.probes)  #number of all non-CGI/non-shore/probes


  ## identification of DMPs in each categolies (TSS, CpG Island, CpG Shore, others)
  TSS.probes <- target_IDs[target_IDs %in% all.TSS.probes]  #identification of TSS DMPs
  CGI.probes <- target_IDs[target_IDs %in% all.CGI.probes]	#identification of CGI DMPs within non-TSS probes
  shore.probes <- target_IDs[target_IDs %in% all.shore.probes]	#identification of shore DMPs within non-TSS probes
  shelf.probes <- target_IDs[target_IDs %in% all.shelf.probes] #identification of shelf DMPs within non-TSS probes
  other.probes <- target_IDs[target_IDs %in% all.other.probes]  #identification of non-TSS, non-CGI, non-shore probes

  ## number of DMPs in each categilies
  nTSS.probes <- length(TSS.probes) #number of CGI DMPs
  nCGI.probes <- length(CGI.probes)	#number of CGI DMPs
  nshore.probes <- length(shore.probes)	#number of shore DMPs
  nshelf.probes <- length(shelf.probes) #number of shore DMPs
  nnonCGI.probes <- length(other.probes)	#number of non-CGI/non-shore DMPs

  ##ramdom sampling from each categolies
  run.TSS.probes <- all.TSS.probes[floor(runif(nTSS.probes, 1, nall.TSS.probes))] #TSS probes
  run.CGI.probes <- all.CGI.probes[floor(runif(nCGI.probes, 1 ,nall.CGI.probes))]	#CGI probes
  run.shore.probes <- all.shore.probes[floor(runif(nshore.probes, 1 ,nall.shore.probes))]	#shore probes
  run.shelf.probes <- all.shore.probes[floor(runif(nshelf.probes, 1 ,nall.shelf.probes))] #shelf probes
  run.other.probes <- all.other.probes[floor(runif(nnonCGI.probes, 1 ,nall.other.probes))] 	#other probes
  run.probes <- c(run.TSS.probes, run.CGI.probes, run.shore.probes, run.shelf.probes, run.other.probes)
  return (run.probes)
}

