#' Utilities
#' 
#' reload of functions
#' 
#' @name utilities
#' 
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @importFrom minfi getAnnotationObject
#' @importFrom methylumi betas methylated unmethylated pval.detect<-
#' 
#' @usage NA
#' @export
#' 

getAnnotationObject <- minfi::getAnnotationObject
betas <- methylumi::betas
methylated <- methylumi::methylated
unmethylated <- methylumi::unmethylated
"pval.detect<-" <- methylumi::"pval.detect<-"