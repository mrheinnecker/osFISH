#!R

library(parallel)
library(getopt)
library(Biostrings)
library(tidyverse)
library(DECIPHER)


spec = matrix(c(
#  'verbose', 'v', 2, "integer",
#  'help' , 'h', 0, "logical",
  'file_target' , 't', 1, "character",
  'file_reference' , 'r', 1, "character",
  'temp_min' , 'a', 1, "double",
  'temp_max' , 'b', 1, "double",
  'length' , 'l', 1, "double",
  'length_range' , 'm', 1, "double"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

source("/g/schwab/Marco/repos/sabeRprobes/R/fncts.R")

# opt <- tibble(
#   file_target="/g/schwab/Marco/projects/osFISH/test2/Prorocentrum_micans_seq_all.fasta",
#   file_reference="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta",
#   temp_min=55,
#   temp_max=75,
#   length=40,
#   length_range=3
# )

t0 <- Sys.time()

message("launching probe design")

probe_length <- as.numeric(opt$length)+as.numeric(opt$length_range)
rel_seq <- readDNAStringSet(opt$file_target)
ref_seq <- readDNAStringSet(opt$file_reference)

cons_seq <- make_consensus_sequence(rel_seq)

temp_center <- mean(c(as.numeric(opt$temp_min), as.numeric(opt$temp_max)))


seq_to_remove <- names(rel_seq)

# Filter out the sequence
filtered_seq <- ref_seq[!(names(ref_seq) %in% seq_to_remove)]

if(length(filtered_seq)==length(ref_seq)){
  message("target sequence not found in reference file")
} else {
  message("removed target sequence from reference file by name match")
}

char_seq <- cons_seq$base %>% paste(collapse="")

all_possible_probes <- mine_probes_from_sequence(char_seq, probe_length, opt$length, opt$length_range) 


test <- DNAStringSet(all_possible_probes%>% pull(2)) 

test <- all_possible_probes %>%
  mutate(
    tm=get_oligo_tm(seq)
  ) %>%
  
  rowwise() %>%
  mutate(len=nchar(seq)) %>%
  mutate(tm_old=64.9 + 41*(str_count(seq, "C|G") - 16.4) / len,
         
         diff_to_center=abs(tm-temp_center))

get_oligo_tm <- function(dna_string){
  
  MeltDNA(dna_string, temps=seq(1, 100)) %>% apply(., 2, function(col){
  
    which.max(col)
  
  })
  
}


#%>%
  mutate(
    meltDNA=
  )





all_possible_probes <- lapply(seq(1, length(rel_seq)), function(nSEQ){
  print(paste(nSEQ, "/",length(rel_seq)))
  char_seq <- as.character(rel_seq[[nSEQ]])
  
  lapply(seq(1, nchar(char_seq)-probe_length), function(START){
    #cat(START, "\n")
    tibble(start=START,
          seq=lapply(seq(as.numeric(opt$length)-as.numeric(opt$length_range),
                         as.numeric(opt$length)+as.numeric(opt$length_range)),
                     function(x){
                       str_sub(char_seq, START, START+x-1)
                     }
                        
            
          ) %>% unlist()
           
    ) %>%
      return()
    
  }) %>% bind_rows() %>%
    mutate(fasta_entry=nSEQ) %>%
    bind_rows() %>%
    return()
  
}) %>%
  bind_rows() %>%
  arrange(fasta_entry, start) %>%
  group_by(seq) %>%
  summarize(
    start=start[1],
    fasta_entry=fasta_entry[1],
    all_fasta_entry=paste(fasta_entry, collapse = ", ")
  ) %>%
  
  rowwise() %>%
  mutate(len=nchar(seq)) %>%
  mutate(tm=64.9 + 41*(str_count(seq, "C|G") - 16.4) / len,
         
         diff_to_center=abs(tm-temp_center)) %>%
  ## filter within each fasta entry (sep. sequence)
  group_by(fasta_entry, start) %>%
  filter(diff_to_center==min(diff_to_center)) %>%
  filter(len==max(len)) %>%
  ungroup() %>%
  mutate(id=paste0("p", seq(1, nrow(.))))
           

remapped <- mclapply(seq(1, nrow(all_possible_probes)), function(n){
#remapped <- lapply(seq(1, 10), function(n){
  
    
  pattern <- all_possible_probes[[n,"seq"]]
  
  print(paste(n, "/", nrow(all_possible_probes)))
  
  target_hits <- rel_seq[which(elementNROWS(vmatchPattern(pattern, rel_seq, max.mismatch = 0))>0)] %>% names()
  target_hits_mm1 <- rel_seq[which(elementNROWS(vmatchPattern(pattern, rel_seq, max.mismatch = 1))>0)] %>% names()
  matches_0 <- filtered_seq[which(elementNROWS(vmatchPattern(pattern, filtered_seq, max.mismatch = 0))>0)] %>% names()
  matches_1 <- filtered_seq[which(elementNROWS(vmatchPattern(pattern, filtered_seq, max.mismatch = 1))>0)] %>% names()
  matches_2 <- filtered_seq[which(elementNROWS(vmatchPattern(pattern, filtered_seq, max.mismatch = 2))>0)] %>% names()
  #matches_2 <- ref_seq[which(elementNROWS(vmatchPattern(pattern, ref_seq, max.mismatch = 2))>0)] %>% names()
  return(
    tibble(id=all_possible_probes[[n,"id"]],
           th_mm0=length(target_hits),
           th_mm1=length(target_hits_mm1),
           th_mm0_perc=length(target_hits)/length(rel_seq),
           th_mm1_perc=length(target_hits_mm1)/length(rel_seq),
           rh_mm0=length(matches_0),
           rh_mm1=length(matches_1),
           rh_mm2=length(matches_2),
           rh_mm0_perc=length(matches_0)/length(filtered_seq),
           rh_mm1_perc=length(matches_1)/length(filtered_seq),
           rh_mm2_perc=length(matches_2)/length(filtered_seq),
           )
  )
}, mc.cores=min(10, detectCores()-1)) %>%
  bind_rows() %>%
  left_join(all_possible_probes %>% select(-diff_to_center))


write_tsv(remapped, file=file.path(str_replace(opt$file_target, "_seq_all.fasta", "_probes.tsv")))

t1 <- Sys.time()

print(paste("Finished in:",t1-t0))
