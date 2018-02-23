writeSplitSeq <- function(seqs, split_num, output_file){
  n_target_seq <- length(seqs)
  ref_point <- 0
  all_filenames <- NULL
  while(ref_point < n_target_seq){
    seqsplit_str <- 1 + ref_point	#start number of sequences
    seqsplit_end <-seqsplit_str + split_num -1	#end number of sequences
    if(seqsplit_end  > n_target_seq){ #for the last incomplete number set
      seqsplit_end <-n_target_seq
    }
    split_seqs <- seqs[seqsplit_str:seqsplit_end]
    filename <- paste(tempDir,"/",output_file,"_",seqsplit_str, "-", seqsplit_end,".fa.tmp" , sep="")
    all_filenames <- c(all_filenames, filename)
    writeXStringSet(split_seqs, file=filename, format="fasta")
    ref_point <- ref_point+split_num  #range sliding
  }
  return(all_filenames) #to read the written files return the file names
}
