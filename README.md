InfiniumDiffMetMotR
===================
Version: 1.0

Description: This is a R package to analyze transcription factor binding motif enrichment for differentially methylated regions.  

Last Update: 2018-3-01

Updated by : takahiro.suzuki.aa@riken.jp

Example
-------
#### 1. motif database construction  
```
library("MotifDb")
targetDB <- "JASPAR_CORE"
targetORG <- c("Hsapiens", "Mmusculus")
targetTF <- "SPI1"
motifDB <- query(MotifDb, targetDB)        #extraction of motif list of "JASPER_CORE"
motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))        #extraction of motifs of "Hsapiens" and "Mmusclus"
motifDB <- query(motifDB,targetTF)       #Extraction of motifs for target TF(s)
motifDBList <- as.list(motifDB)
```

#### 2. Screening of enriched motifs  
```
library(InfiniumDiffMetMotR)
MotScr(infile="sel_processed_Mval.txt", motifDBList = motifDBList, cutoff=2, outname="screening_result", ControlColnum=1, TreatmentColnum=2, MethylDemethyl="Demethyl", version = "EPIC")
```
infile: M-value matrix of infinium methylation array
motifDBList: list. Each factor is PWM as fllowing format;
  1 2 3         4         5         6          7          8          9
A 0 0 0 0.1189189 0.1027027 0.2972973 0.28648649 0.10270270 0.04864865
C 0 1 1 0.3837838 0.3081081 0.2378378 0.16216216 0.08648649 0.42162162
G 1 0 0 0.2486486 0.3297297 0.3621622 0.49189189 0.74054054 0.42702703
T 0 0 0 0.2486486 0.2594595 0.1027027 0.05945946 0.07027027 0.10270270
cutoff: cutoff velue of delta M.
outname: name of output file
ControlColnum: A column of control data, such as unddiferentiated.
TreatmentColnum: A column of treatment data such as differentiated.
MethylDemethyl: "Methyl" (or "M") or "dimethyl (or"D").analysis target of differentially methylated probes.
 version: "450" or "EPIC" (or "850"). version of methylation array.

