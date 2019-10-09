#' GetGEO idat file & M-value conversion
#' 
#' Downloads idat files based given GEO ID, normalzies, and computes M-values.
#' 
#' @param GEOid GEO ID
#' @param version type of array version for analysis. 
#' @param sampleNames a vector of sampeNamse. If "FALSE", sample namses use GEO sample titles.
#' 
#' @importFrom GEOquery getGEO Meta getGEOSuppFiles gunzip
#' @importFrom utils untar
#' 
#' @return Normalized M-value text files List.
#' @keywords GEO idat Normalization M-value
#' @export


getIdat2M <- function(GEOid, version, sampleNames=FALSE){
    original_dir <- getwd()
    getGEOSuppFiles(GEOid)    #get supplimentary files

    setwd(paste0(GEOid, "/"))    #move to the supplimenrary file directory

    untar(tarfile = paste0(GEOid, "_RAW.tar"))

    ##Uncompress the idat files
    idat_files <- list.files(pattern = "idat.gz")
    
    for(i in 1:length(idat_files)){
        gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
    }

    ## idat files are collected to a directory by version
    idat_unzips <- list.files(pattern = "idat$")
    idat_GEOid <- sapply(strsplit(idat_unzips, "_"), function(x){x[1]})
	info_GEOList <- lapply(idat_GEOid, getGEO)
    platform <- sapply(info_GEOList, function(x){Meta(x)$platform_id})
    methyl27_index <- grep("GPL8490", platform)
    methyl450_index <- grep("GPL13534", platform)
    EPIC_index <- grep("GPL21145|GPL23976", platform)

    cat(paste(length(EPIC_index), "HumanMethylationEpic idat files found", sep=" "))
    cat("\n")
    cat(paste(length(methyl450_index), "HumanMethylation450 idat files found", sep=" "))
    cat("\n")
    cat(paste(length(methyl27_index), "HumanMethylation27 idat files found", sep=" "))
    cat("\n")
    
    if(!identical(EPIC_index, integer(0))){
        dir.create("EPIC")
        for(i in idat_unzips[EPIC_index]){
            file.rename(i, paste0("EPIC/", i))
        }
    }
    if(!identical(methyl450_index, integer(0))){
        dir.create("Methyl450")
        for(i in idat_unzips[methyl450_index]){
            file.rename(i, paste0("Methyl450/", i))
        }
    }
    if(!identical(methyl27_index, integer(0))){
        dir.create("Methyl27")
        for(i in idat_unzips[methyl27_index]){
            file.rename(i, paste0("Methyl27/", i))
        }
    }
    if(identical(EPIC_index, integer(0)) && identical(methyl450_index, integer(0)) && identical(methyl27_index, integer(0))){
        stop("Methylation array not found")
    }

    ## Normalization
    selDataMatrixList <- as.list(NULL)
    for (j in version){
        v.dir <- switch(j,
            "EPIC" = "EPIC",
            "850" = "EPIC",
            "450" = "Methyl450",
            "27" = "Methyl27",
            stop("Array varsions not correct")
            )
        setwd(paste0(v.dir, "/"))
        cat(paste("Normalizing", j, "Arrays...", sep=" "))
        cat("\n")
        selDataMatrix <- lumiMethyNorm(idatpath=getwd(), inputtype = "idat")    #Normalization and M-value computation
        if(length(version) != 1){
            selDataMatrixList <- c(selDataMatrix, v.dir=list(selDataMatrix))
        }
        setwd("../")
    }
    setwd(original_dir)

    ## Return object depends on indicated version length
    if (length(version) == 1){
        return_object <- selDataMatrix
    }else{
        return_object <- selDataMatrixList
    }
    return(return_object)
}


