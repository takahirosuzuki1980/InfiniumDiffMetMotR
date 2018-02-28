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
`motifDB <- query(MotifDb, targetDB)`        #extraction of motif list of "JASPER_CORE"  
`motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))`        #extraction of motifs of "Hsapiens" and "Mmusclus"  
`motifDB <- query(motifDB,targetTF)`       #Extraction of motifs for target TF(s)  
`motifDBList <- as.list(motifDB)`  
##discription of motifs  
`motGene <- values(motifDB)[,4]`  
`motID <- values(motifDB)[,2]`  
`motSource <- values(motifDB)[,3]`  
`motOrg <- values(motifDB)[,9]`  

#### 2. Identification of differentially methylated regions
`infile <- sel_processed_Mval.txt`
`outname <-iPS-HPC_SPI1`
`ControlColnum <- 1`
`TreatmentColnum <- 11`
`MethylDemethyl <- "Demethyl"`

`selDataMatrix <- read.table (infile)`
#extraction demethylated probes
`if((MethylDemethyl == "Demethyl") ||( MethylDemethyl == "D")) {`
`	diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) >=2)`
`}else if ((MethylDemethyl == "Methyl" )|| (MethylDemethyl == "M")){`
`	diff_table <- which((selDataMatrix[,ControlColnum]-selDataMatrix[,TreatmentColnum]) <=-2)`
`}`

`DMP_IDs <- rownames(selDataMatrix )[diff_table]`
`nDMP_IDs <- length(DMP_IDs)`
`allProbe_IDs <- rownames(selDataMatrix)`
