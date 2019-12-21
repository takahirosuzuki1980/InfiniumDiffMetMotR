#' Write sequences to multiple output files
#' 
#' To reduce the memory occupancy, the DNAStringSet sequences are output to files, splitted into indicate numbers.
#' 
#' @param seqs DNAStringSet sequences
#' @param split_num nuber of sequences per file
#' @param tempDir directory to be stored the sequences
#' @param output_file file ID for output files
#' 
#' @importFrom Biostrings writeXStringSet
#' 
#' @return output filenames
#' 
#' @keywords sequences DNAStringSet split
#' @export

writeSplitSeq <- function(seqs, split_num=2500, tempDir="Transient.out", output_file="transient.out"){
  n_target_seq <- length(seqs)
  ref_point <- 0
  all_filenames <- NULL
  while(ref_point < n_target_seq){
    seqsplit_str <- 1 + ref_point    #start number of sequences
    seqsplit_end <-seqsplit_str + split_num - 1    #end number of sequences
    if(seqsplit_end  > n_target_seq){    #for the last incomplete number set
      seqsplit_end <-n_target_seq
    }
    split_seqs <- seqs[seqsplit_str:seqsplit_end]
    out_name <- paste(tempDir,"/",output_file,"_",seqsplit_str, "-", seqsplit_end,".fa.tmp" , sep="")
    all_filenames <- c(all_filenames, out_name)
    writeXStringSet(split_seqs, filepath=out_name, format="fasta")
    ref_point <- ref_point+split_num    #range sliding
  }
  return(all_filenames)    #to read the written files return the file names
}
