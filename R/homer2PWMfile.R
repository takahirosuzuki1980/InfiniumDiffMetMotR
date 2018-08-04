homer2PWMfile <- function(file){
	table <- read.table(file=file, fill=T, row.names=NULL, header=F)
	table.vector <- cbind(as.vector(table[[1]]), as.vector(table[[2]]), as.vector(table[[3]]), as.vector(table[[4]]))

	IDs <- grep(">",table.vector[,1])
	names <- table.vector[IDs,2]
	listPWM <- list()

	for (i in seq(length(IDs))){
    if (i == length(IDs)){ ##for the last motif
        table.vector.end <- nrow(table.vector)
    }else{
        table.vector.end <- (IDs[i+1]-1)
    }
	listPWM[[i]] <- matrix(as.numeric(table.vector[(IDs[i]+1):table.vector.end,1:4]), nrow=4, byrow=T, dimnames=list(c("A","C","G","T")))
	colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
	names(listPWM)[i] <-  names[i]
	}
	return (listPWM)
}