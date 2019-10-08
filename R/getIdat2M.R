#' GetGEO idat file & M-value conversion
#' 
#' Downloads idat files based given GEO ID, normalzies, and computes M-values.
#' 
#' @param GEOid GEO ID
#' @param sampleNames a vector of sampeNamse. If "FALSE", sample namses use GEO sample titles.
#' 
#' @importFrom GEOquery getGEOSuppFiles gunzip
#' @importFrom utils untar
#' 
#' @return Normalized M-value text files
#' @keywords GEO idat Normalization M-value
#' @export

getIdat2M <- function(GEOid, sampleNames=FALSE){
    getGEOSuppFiles(GEOid)    #get supplimentary files

    setwd(paste0(GEOid, "/"))    #move to the supplimenrary file directory

    untar(tarfile = paste0(GEOid, "_RAW.tar"))

    ##Uncompress the idat files
    idat_files <- list.files(pattern = "idat.gz")
    
    for(i in 1:length(idat_files)){
        gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
    }

    lumiMethyNorm(idatpath=getwd(), inputtype = "idat")    #Normalization and M-value computation
    setwd("../")
}


