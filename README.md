InfiniumDiffMetMotR
===================
Version: 1.3

Description: This is a R package to analyze transcription factor binding motif enrichment at differentially methylated regions for Infinium Methylation BeadChip (Illumina).  

Last Update: 2018-06-12  

Depends: R (>= 2.10), Biobase (>= 2.5.5)  

Author: Takahiro Suzuki  

Updated by: takahiro.suzuki.aa@riken.jp

Install
-------
Install of dependent Bioconductor packages
```
source("https://bioconductor.org/biocLite.R")
biocLite("FDb.InfiniumMethylation.hg19")
biocLite("lumi")
biocLite("Biostrings")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
```
Install of InfiniumDiffMetMotR from github
```
install.packages("devtools")
require(devtools)
install_github("takahirosuzuki1980/InfiniumDiffMetMotR")
```

Example
-------
#### 1. Load of InfiniumDiffMetMotR
```
library("InfiniumDiffMetMotR")
```

#### 2. Normalization  
```
library(InfiniumDiffMetMotR)  
lumiMethyNorm(fileName = "TableControl.txt", sample_names = sample_names)
```
input unnourmalized-no background correction full infinium methylation array data.  
The header of the data should includes the following items as the following formats;  
  
**TargetID    Sample1.Signal_A    Sample1.Signal_B    Sample1.Detection Pval**  
  
TargetID: Illumina ID  
Sample1.Signal_A: Signal of Unmethylated probe of sample 1  
Sample1.Signal_B: Signal of Methylated probe of sample 1  
Sample1.Detection Pval: detection P value of sample 1  
  
sample_names: vector of sample names  
  
lumiMethyNorm generates **Process_Result folder**: reports of the normalization, **processed_Mval.txt**: a matrix of M-value, and **sel_processed_Mval.txt**: a matrix of M-value which does not inculude low detection p-value probes. (default cut off is 0.01)  
  
#### 3. motif database construction  
**Example 1: JASPER_CORE, Hsapiencs and Mmusclus**
```
library("MotifDb")
targetDB <- "JASPAR_CORE"
targetORG <- c("Hsapiens", "Mmusculus")
motifDB <- query(MotifDb, targetDB)        #extraction of motif list of "JASPER_CORE"
motifDB <- c(query(motifDB,targetORG[1]),query(motifDB,targetORG[2]))        #extraction of motifs of "Hsapiens" and "Mmusclus"
```
If you want analyze a specific motif select a motif. (ex. SPI1)
```
targetTF <- "SPI1"
motifDB <- query(motifDB,targetTF)       #Extraction of motifs for target TF(s)
```
Finally, convert the motif list to list format.
```
motifDBList <- as.list(motifDB)
```
**Example 2: IMAGE motif database**
```
motifDBList <- IMAGE_PWMlist
```
If you want to extract the motif of interest (e.g. SPI1):
```
motif_names <- "SPI1"
motifDBList <- motifDBList[grep(motif_names,names(motifDBList))]
```
Motif list should be a list of PWMs of following format;  
  
   1     2    3    4    5    6    7    8    9  
A    0    0    0    0.1189189    0.1027027    0.2972973    0.28648649    0.10270270    0.04864865  
C    0    1    1    0.3837838    0.3081081    0.2378378    0.16216216    0.08648649    0.42162162  
G    1    0    0    0.2486486    0.3297297    0.3621622    0.49189189    0.74054054    0.42702703  
T    0    0    0    0.2486486    0.2594595    0.1027027    0.05945946    0.07027027    0.10270270  
  
 
#### 4. Screening of enriched motifs  
```
MotScr(infile="sel_processed_Mval.txt", motifDBList = motifDBList, cutoff = 2, p.cutoff = 0.001, outname="screening_result", ControlColnum=c(1,2), TreatmentColnum=c(3,4), MethylDemethyl="Demethyl", version = "850")
```
infile: M-value matrix of infinium methylation array  
motifDBList: list format of PWM  
cutoff: cutoff velue of delta M  
outname: name of output file  
ControlColnum: A column of control data, such as unddiferentiated.  
TreatmentColnum: A column of treatment data such as differentiated.  
MethylDemethyl: "Methyl" (or "M") or "dimethyl (or"D").analysis target of differentially methylated probes.  
version: "450" or "EPIC" (or "850"). version of methylation array.  
  
If you perform a comparison of multiple samples, identification of differentially methylated probes uses both Welch's t-test and M-value difference. For single sample comparison, it uses only M-value difference (delta M).

