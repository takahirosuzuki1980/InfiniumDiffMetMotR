lumiMethyNorm <- function(fileName = "TableControl.txt", sample_names = sample_names){
	################ ライブラリーを読み込み #########################
	library(lumi)
	library(annotate)
	library(FDb.InfiniumMethylation.hg19)
	############### データの読み込み #############################
	data.lumiMethy <- lumiMethyR(fileName)
	addAnnotationInfo(data.lumiMethy, lib = 'FDb.InfiniumMethylation.hg19', annotationColumn=c('COLOR_CHANNEL', 'CHROMOSOME', 'POSITION'))
	sampleNames(data.lumiMethy) <- sample_names	#convert the sampleID to sample name
	############### 必要なら###################################
	dir.create("Process_Result")	#make a directory to store processing results
	##PCA plot
	pdf("Process_Result/pca.pdf")
	plotSampleRelation(data.lumiMethy, method='mds', cv.Th=0)
	dev.off()
	################# Color balance adjustment between two color channels ###########################
	lumiMethy.c.adj <- lumiMethyC(data.lumiMethy)
	## Check color balance after color balance adjustment
	pdf("Process_Result/col.adj.pdf")
	par(mfrow=c(2,2))
	plotColorBias1D(data.lumiMethy, channel='sum')	#density plot of raw
	boxplotColorBias(data.lumiMethy, channel='sum')	#boxplot of raw
	plotColorBias1D(lumiMethy.c.adj, channel='sum')	#density plot of adjusted
	boxplotColorBias(lumiMethy.c.adj, channel='sum')	#box plot of adjusted
	dev.off()
	############################### Normalization###############################
	## Perform SSN normalization based on color balance adjusted data
	## perform quantile normalization based on color balance adjusted data
	lumiMethy.c.q <- lumiMethyN(lumiMethy.c.adj, method='quantile')
	## plot the density of M-values before and after quantile normalization
	pdf("Process_Result/normalize.pdf")
	par(mfrow=c(2,2))
	density(data.lumiMethy, main="Density plot of M-value after quantile normalization")
	plotSampleRelation(data.lumiMethy, method='cluster', cv.Th=0)
	density(lumiMethy.c.q, main="Density plot of M-value after quantile normalization")
	plotSampleRelation(lumiMethy.c.q, method='cluster', cv.Th=0)
	dev.off()
	###################### output the normlized M-value as a Tab-separated text file ###########################
	dataMatrix <- exprs(lumiMethy.c.q)
	write.exprs(lumiMethy.c.q, file='processed_Mval.txt')
	#################### Remove the probes which is low detection p-value among all sample.##########################
	presentCount <- detectionCall(data.lumiMethy, Th = 0.01)	#Sample number of "Present" of a probe.
	selDataMatrix <- dataMatrix[presentCount > 0,]
	write.table(selDataMatrix, file="sel_processed_Mval.txt", sep="\t", quote=F)
}
