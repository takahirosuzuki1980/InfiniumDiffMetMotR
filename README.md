InfiniumDiffMetMotR
===================
This is a R package to analyze transcription factor binding motif enrichment for differentially methylated regions.  

Example
-------
#### 1. motif database construction  
`library("MotifDb")`
`targetDB <- "JASPAR_CORE"`  
`targetORG <- c("Hsapiens", "Mmusculus")`  
`targetTF <- "SPI1"`  
`motifDB <- query(MotifDb, targetDB)`        #extract motif list of "JASPER_CORE"  
`motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))`        #extraction of motifs of "Hsapiens" and "Mmusclus"  
`motifDB <- query(motifDB,targetTF)`       #Extraction of motifs for target TF(s)  
`motifDBList <- as.list(motifDB)`  
##discription of motifs  
`motGene <- values(motifDB)[,4]`  
`motID <- values(motifDB)[,2]`  
`motSource <- values(motifDB)[,3]`  
`motOrg <- values(motifDB)[,9]`  
