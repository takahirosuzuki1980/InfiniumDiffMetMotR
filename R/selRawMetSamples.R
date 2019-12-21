#' data selection for raw Infinium signal data
#' 
#' This function extracts necessary data of targets from the infinium array raw data table.
#' 
#' @param TableControl  unnormalized no background corrction raw methylation array data
#' @param sampleIDs a vector of sample IDs to be selected
#' 
#' @importFrom dplyr %>% select one_of
#' @importFrom utils read.table write.table
#' 
#' @return selected unnormalized no background corrction raw methylation array data
#' @keywords Infinium, sampleID
#' @export
#' 
selRawMetSamples <- function(TableControl = "TableControl.txt", sampleIDs){
    raw_table <- read.table(TableControl, sep="\t", header=T, stringsAsFactors=F)    # ファイルの読み込む
    IDs <-  sampleIDs    # サンプルIDの読み込み

    ## ID conversion (それぞれのIDにSignal_A, Signal_B, Detection.Pvalを付け加える)
    ## Rでは列名の最初に数字が使えないのでIDが数字の場合データを読み込んだ時に自動的にXがはいる。その場合はXも最初に付け加えている、
    ID_header <- "TargetID"
    for (i in IDs){
        if(grepl("^[0-9]", i)){    # IDが数字の場合
            i <- paste0("X", i)    # IDにXをつける
        }
        ID_header <- c(ID_header, paste0(i, ".Signal_A"))
        ID_header <- c(ID_header, paste0(i, ".Signal_B"))
        ID_header <- c(ID_header, paste0(i, ".Detection.Pval"))
    }
    
    sel_raw_table <- raw_table %>% select(one_of(ID_header))    # IDに一致する列を抽出
    write.table(sel_raw_table, file="sel_TableControl.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)    #データの書き出し
}