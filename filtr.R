library(getopt)
library(Biostrings)
library(tidyverse)


spec = matrix(c(
  #  'verbose', 'v', 2, "integer",
  #  'help' , 'h', 0, "logical",
  'file_reference' , 'r', 1, "character",
  'outdir' , 'a', 1, "character",
  'target' , 'b', 1, "character",
  
), byrow=TRUE, ncol=4)


opt = getopt(spec)

opt <- tibble(
  target="/g/schwab/Marco/projects/osFISH/RCC_cultures_ordered.txt",
  file_reference="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta",
  outdir="/g/schwab/Marco/projects/osFISH/test/", 
  store_separate=F
  # temp_max=75, 
  # length=40,
  # length_range=3
)


extract_by_str_match <- function(ref_seq, target, exact_name=FALSE){
  
  ooi <-  tolower(target)
  seq_ooi <- ref_seq[which(str_detect(tolower(names(ref_seq)), ooi))]
  
  return(seq_ooi)
}


#align_multiple_target_sequences <- function(all_seq_matches)






ref_seq <- readDNAStringSet(opt$file_reference)

if(file.exists(opt$target)){
  
  target <- read_tsv(opt$target, col_names = FALSE) %>% pull(1)
  
} else {
  
  target <- opt$target
  
}

#TARGET <- "Heterocapsa_rotundata"

combined_outfiles <- lapply(target, function(TARGET){
  
  message(TARGET)
  all_seq_matches <- extract_by_str_match(ref_seq, TARGET)
  
  message("  ",length(all_seq_matches), " sequences found")
  
  outdir <- opt$outdir
  
  if(length(all_seq_matches)==0){
    warning("no matches found")
    all_outfiles <- tibble()
  } else if(opt$store_separate){
    message("storing all", length(all_seq_matches), "sequences in separate files at:", outdir)
    all_outfiles <- lapply(seq(1,length(all_seq_matches)), function(n){
      
      outfile <- file.path(outdir, paste0(TARGET, "_seq", n, ".fasta"))
      
      writeXStringSet(all_seq_matches[n], 
                      file = outfile)
      
      return(tibble(file=outfile))
      
    }) %>%
      bind_rows()
    
  } else {
    outfile <- file.path(outdir, paste0(TARGET, "_seq_all.fasta"))
    writeXStringSet(all_seq_matches, 
                    file = outfile)
    all_outfiles <- tibble(file=outfile)
  }
  
    return(all_outfiles)

}) %>%
  bind_rows()


write_csv(combined_outfiles,
          file=file.path(outdir, "output_paths.tsv"))









