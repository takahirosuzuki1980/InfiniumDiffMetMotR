#' generation of GEO sbmission files
#' 
#' This scripts generate Matrix_signal_intensities.txt and M-value based Matrtix_processed.txt for GEO submission.
#' 
#' @param idat_path path to idat files
#' @param sample_names a vactor of sample names
#' 
#' @importFrom minfi read.metharray.exp preprocessIllumina detectionP getMeth getUnmeth
#' @importFrom dplyr %>% mutate
#' 
#' @return a matrix containg number of PWM at demethylated regions, number of PWM at random regions, TF, and adjust p-value
#' 
#' @keywords count of PWM, fold-cange, p-value
#' @export
#' 

idat2GEO <- function(idat_path, sample_names = NULL){
    #Unnormalized signal data (Matrix_signal_intensities.txt)
    RGset <- read.metharray.exp(base = idat_path)
    NNdata <- preprocessIllumina(RGset, bg.correct = FALSE, normalize = "no")

    detection.p <- detectionP(RGset)
    colnames(detection.p) <- paste0(colnames(detection.p), ".Detection Pval")
    detection.p  %>% 
        mutate(ID_REF = rownames(detection.p)) -> detection.p

    m_signal <- getMeth(NNdata)
    colnames(m_signal) <- paste0(colnames(m_signal), ".Methylated_signal")
    m_signal  %>% 
        mutate(ID_REF = rownames(m_signal)) -> m_signal

    u_signal <- getUnmeth(NNdata)
    colnames(u_signal) <- paste0(colnames(u_signal), ".Unmethylated_signal")
    u_signal %>%
        mutate(ID_REF = rownames(u_signal)) -> u_signal

    #ID_REF, methylated signal, unmethylated signal, detection p-valueを統合
    #sample_namesがない場合を考慮
    #書き出し


    #Normalized M-value data (Matrix_processed.txt)

}