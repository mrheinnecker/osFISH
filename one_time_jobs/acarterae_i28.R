#!R/4.4.0
## align imager28 to Amphidinium carterae rRNA and or genome


#BiocManager::install("pwalign")
.libPaths("/home/rheinnec/R/x86_64-pc-linux-gnu-library/4.4")
library(tidyverse)
library(Biostrings)
# helper: count matches/mismatches/gaps from aligned strings
.aln_stats <- function(p_aln, s_aln) {
  p <- as.character(p_aln)
  s <- as.character(s_aln)
  
  p_ch <- str_split(p, "", simplify = TRUE)
  s_ch <- str_split(s, "", simplify = TRUE)
  
  stopifnot(ncol(p_ch) == ncol(s_ch))
  
  gap_p <- p_ch == "-"
  gap_s <- s_ch == "-"
  gaps  <- sum(gap_p | gap_s)
  
  comparable <- !(gap_p | gap_s)
  matches    <- sum(comparable & (p_ch == s_ch))
  mismatches <- sum(comparable & (p_ch != s_ch))
  aligned_len <- ncol(p_ch)
  
  tibble(
    aligned_width = aligned_len,
    matches = matches,
    mismatches = mismatches,
    gaps = gaps,
    identity = ifelse(matches + mismatches > 0, matches / (matches + mismatches), NA_real_)
  )
}

# main: list(alignment objects) -> tibble
alignments_to_tibble <- function(aln_list, keep_strings = TRUE) {
  imap_dfr(aln_list, function(aln, idx) {
    # These accessors work for PairwiseAlignments objects
    pat_start <- start(pattern(aln))
    pat_width <- width(pattern(aln))
    pat_end   <- pat_start + pat_width - 1
    
    sub_start <- start(subject(aln))
    sub_width <- width(subject(aln))
    sub_end   <- sub_start + sub_width - 1
    
    sc <- score(aln)
    
    p_aln <- alignedPattern(aln)
    s_aln <- alignedSubject(aln)
    
    stats <- .aln_stats(p_aln, s_aln)
    
    tibble(
      list_index = idx,
      score = as.numeric(sc),
      
      pattern_start = as.integer(pat_start),
      pattern_end   = as.integer(pat_end),
      pattern_width = as.integer(pat_width),
      
      subject_start = as.integer(sub_start),
      subject_end   = as.integer(sub_end),
      subject_width = as.integer(sub_width)
    ) %>%
      bind_cols(stats) %>%
      mutate(
        pattern_seq = if (keep_strings) as.character(pattern(aln)) else NA_character_,
        subject_seq = if (keep_strings) as.character(subject(aln)) else NA_character_,
        aligned_pattern = if (keep_strings) as.character(p_aln) else NA_character_,
        aligned_subject = if (keep_strings) as.character(s_aln) else NA_character_
      )
  })
}


opt <- tibble(
  #file_target="/g/schwab/marco/projects/osFISH/test2/Prorocentrum_micans_seq_all.fasta",
  file_reference="/g/schwab/marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta",
  #fullseq_fasta="/g/schwab/marco/projects/osFISH/JBNFNT01.1.fsa_nt.gz"
  fullseq_fasta="/g/schwab/marco/projects/osFISH/fasta/GSE94355_expressed.fa.gz"
)

sm <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)


full_seq <- readDNAStringSet(opt$fullseq_fasta)


# ref_seq <- readDNAStringSet(opt$file_reference)
# 
# 
# ac_seqs <- ref_seq[which(str_detect(names(ref_seq), "Amphidinium_carterae"))]

i28_seq <- DNAStringSet("TTGTTAAGTTGTGTTAAGTTGTTTTTTTT") %>% reverseComplement()
i25_seq <- DNAStringSet("TTTATTATTGGTATTATTGGTTTTTTTT") %>% reverseComplement()
i27_seq <- DNAStringSet("TTATGATGATGTATGATGATGTTTTTTTT") %>% reverseComplement()
i28_nl <- DNAStringSet("TTGTTAAGTTGTGTTAAGTTGT") %>% reverseComplement()
i25_nl <- DNAStringSet("TTTATTATTGGTATTATTGGT") %>% reverseComplement()
i27_nl <- DNAStringSet("TTATGATGATGTATGATGATGT") %>% reverseComplement()

to_test <- c(i28=i28_seq, 
             i27=i27_seq, 
             i28_nl=i28_nl,
             i27_nl=i27_nl,
             #i25=i25_seq, 
             i28rc=reverseComplement(i28_seq), 
             i27rc=reverseComplement(i27_seq),
             i28_nl_rc=reverseComplement(i28_nl), 
             i27_nl_rc=reverseComplement(i27_nl)
             )
             #i25rc=reverseComplement(i25_seq))



all <- lapply(names(to_test),function(nt){
  
  
  target <- to_test[[nt]]
  print(nt)
  print(target)
  
  alignments <- lapply(seq(1,length(full_seq)), function(ns){
  #  alignments <- lapply(seq(1,5), function(ns){
    print(ns)
    subject <- full_seq[[ns]]
    #print(subject)
    pairwiseAlignment(
      target, 
      subject, 
      type = "local",
      #substitutionMatrix = sm,
      gapOpening = -5,      # less negative => easier to open a gap
      gapExtension = -3     # less negative => easier to extend a gap
      
    )
    
  })
  aln_tbl <- alignments_to_tibble(alignments, keep_strings = TRUE) %>%
    mutate(target=nt)
  write_rds(alignments, file=file.path("/g/schwab/marco/projects/osFISH/alignments", paste0(nt, "_aln_1.rds")))
  write_tsv(aln_tbl, file=file.path("/g/schwab/marco/projects/osFISH/alignments", paste0(nt, "_aln_1.tsv")))
  return(list(alignments, aln_tbl))
  
})



## ribosolmal RNA Acarterae


ref_seq <- readDNAStringSet(opt$file_reference)


ac_seqs <- ref_seq[which(str_detect(names(ref_seq), "Amphidinium_carterae"))]



all_rRNA <- parallel::mclapply(names(to_test),function(nt){
  
  
  target <- to_test[[nt]]
  print(nt)
  print(target)
  
  alignments <- lapply(seq(1,length(ac_seqs)), function(ns){
    #  alignments <- lapply(seq(1,5), function(ns){
    print(ns)
    subject <- ac_seqs[[ns]]
    #print(subject)
    pairwiseAlignment(
      target, 
      subject, 
      type = "local",
      #substitutionMatrix = sm,
      gapOpening = -3,      # less negative => easier to open a gap
      gapExtension = -2     # less negative => easier to extend a gap
      
    )
    
  })
  aln_tbl <- alignments_to_tibble(alignments, keep_strings = TRUE) %>%
    mutate(target=nt)
  write_rds(alignments, file=file.path("/g/schwab/marco/projects/osFISH/alignments", paste0(nt, "_aln_2.rds")))
  write_tsv(aln_tbl, file=file.path("/g/schwab/marco/projects/osFISH/alignments", paste0(nt, "_aln_2.tsv")))
  return(list(alignments, aln_tbl))
  
})

ov_rRNA <- lapply(all_rRNA, dplyr::first)

full_rRNA <- lapply(all_rRNA, last) %>%
  bind_rows()

full_rRNA %>%
  arrange(desc(score)) %>%
  select(1,2,9:18) %>%
  View()




files <- list.files("/g/schwab/marco/projects/osFISH/alignments/", pattern=".tsv", full.names=T)

full <- lapply(files, function(FILE){
  
  read_tsv(FILE) %>%
    mutate(file=basename(FILE))
  
}) %>%
  bind_rows()



ov <- full %>%
  group_by(file) %>%
  arrange(desc(score)) %>%
  mutate(id=c(1:length(file))) %>%
  filter(id %in% 1:3, !str_detect(file, "rc")) %>%
  select(1,2,9:15)

check_i28 <- read_rds("/g/schwab/marco/projects/osFISH/alignments/i28_aln.rds")
check_i27 <- read_rds("/g/schwab/marco/projects/osFISH/alignments/i27_aln.rds")


transcripts <- tibble(
  list_index=ov$list_index,
  name=names(full_seq[ov$list_index])
) %>% unique()


dfa <- readxl::read_xlsx("/g/schwab/marco/projects/osFISH/alignments/GSE94355_data_all_expressed.corr (1).xlsx") %>%
  rownames_to_column("list_index") 
#final <- 
  ov %>% select(1,2,3,5,6,7,8,9,10) %>%
  unique() %>% 
  left_join(transcripts, by="list_index") %>%
  mutate(list_index=as.character(list_index)) %>%
  left_join(
    dfa %>%select(1:5, list_index),
    by=c("list_index")
  ) %>%
  ungroup() %>%
  mutate(
    mean_counts=(c1.counts+c2.counts+c3.counts)/3,
    seq=str_split(file, "_") %>% map_chr(.,1)
  ) %>%
  select(-c1.counts, -c2.counts, -c3.counts, -file, -name, -id, -list_index) %>%
    unique() %>%
    filter(matches>10)










## 
# files <- list.files("/g/schwab/marco/projects/osFISH/alignments/", pattern=".tsv", full.names=T)
# 
# full <- lapply(files, function(FILE){
#   
#   read_tsv(FILE) %>%
#     mutate(file=basename(FILE))
#   
# }) %>%
#   bind_rows()
# 
# 
# 
# full %>%
#   group_by(file) %>%
#   arrange(desc(score)) %>%
#   mutate(id=c(1:length(file))) %>%
#   filter(id %in% 1:5) %>%
#   select(1,2,9:15) %>%
#   View()






# 
# 
# ac_probe <- DNAStringSet("TTTGCTGCACCCTTCCTCAGCACTTGAGCCTAGAGA")%>% reverseComplement()
# 
# alignments <- lapply(ac_seqs, function(subject){
#   print(1)
#   pairwiseAlignment(
#     i28_seq, 
#     subject, 
#     type = "local",
#     #substitutionMatrix = sm,
#     gapOpening = -3,      # less negative => easier to open a gap
#     gapExtension = -2     # less negative => easier to extend a gap
#     
#   )
#   
# })
# 
# alignments_full <- lapply(seq(1,length(full_seq)), function(ns){
#   print(ns)
#   subject <- full_seq[ns]
#   pairwiseAlignment(
#     i28_seq, 
#     subject, 
#     type = "local",
#     #substitutionMatrix = sm,
#     #gapOpening = -3,      # less negative => easier to open a gap
#     #gapExtension = -2     # less negative => easier to extend a gap
#     
#   )
#   
# })
# 
# 
# alignments_full_revcomp <- lapply(seq(1,length(full_seq)), function(ns){
#   print(ns)
#   subject <- full_seq[ns]
#   pairwiseAlignment(
#     reverseComplement(i28_seq), 
#     subject, 
#     type = "local",
#     #substitutionMatrix = sm,
#     #gapOpening = -3,      # less negative => easier to open a gap
#     #gapExtension = -2     # less negative => easier to extend a gap
#     
#   )
#   
# })
# 
# 
# 
# 
# 

# 
# # helper: count matches/mismatches/gaps from aligned strings
# .aln_stats <- function(p_aln, s_aln) {
#   p <- as.character(p_aln)
#   s <- as.character(s_aln)
#   
#   p_ch <- str_split(p, "", simplify = TRUE)
#   s_ch <- str_split(s, "", simplify = TRUE)
#   
#   stopifnot(ncol(p_ch) == ncol(s_ch))
#   
#   gap_p <- p_ch == "-"
#   gap_s <- s_ch == "-"
#   gaps  <- sum(gap_p | gap_s)
#   
#   comparable <- !(gap_p | gap_s)
#   matches    <- sum(comparable & (p_ch == s_ch))
#   mismatches <- sum(comparable & (p_ch != s_ch))
#   aligned_len <- ncol(p_ch)
#   
#   tibble(
#     aligned_width = aligned_len,
#     matches = matches,
#     mismatches = mismatches,
#     gaps = gaps,
#     identity = ifelse(matches + mismatches > 0, matches / (matches + mismatches), NA_real_)
#   )
# }
# 
# # main: list(alignment objects) -> tibble
# alignments_to_tibble <- function(aln_list, keep_strings = TRUE) {
#   imap_dfr(aln_list, function(aln, idx) {
#     # These accessors work for PairwiseAlignments objects
#     pat_start <- start(pattern(aln))
#     pat_width <- width(pattern(aln))
#     pat_end   <- pat_start + pat_width - 1
#     
#     sub_start <- start(subject(aln))
#     sub_width <- width(subject(aln))
#     sub_end   <- sub_start + sub_width - 1
#     
#     sc <- score(aln)
#     
#     p_aln <- alignedPattern(aln)
#     s_aln <- alignedSubject(aln)
#     
#     stats <- .aln_stats(p_aln, s_aln)
#     
#     tibble(
#       list_index = idx,
#       score = as.numeric(sc),
#       
#       pattern_start = as.integer(pat_start),
#       pattern_end   = as.integer(pat_end),
#       pattern_width = as.integer(pat_width),
#       
#       subject_start = as.integer(sub_start),
#       subject_end   = as.integer(sub_end),
#       subject_width = as.integer(sub_width)
#     ) %>%
#       bind_cols(stats) %>%
#       mutate(
#         pattern_seq = if (keep_strings) as.character(pattern(aln)) else NA_character_,
#         subject_seq = if (keep_strings) as.character(subject(aln)) else NA_character_,
#         aligned_pattern = if (keep_strings) as.character(p_aln) else NA_character_,
#         aligned_subject = if (keep_strings) as.character(s_aln) else NA_character_
#       )
#   })
# }

# # ---- usage ----
# # aln_list <- your list of local alignments (like the one you printed)
# aln_tbl <- alignments_to_tibble(alignments_full, keep_strings = TRUE)
# 
# alignments_full_revcomp
# 
# aln_revcomp <- alignments_to_tibble(alignments_full_revcomp, keep_strings = TRUE)
# 
# 
# # show "best ones" (tune ordering however you like)
# best <- aln_tbl %>%
#   arrange(desc(score), desc(identity), desc(pattern_width), subject_start) %>%
#   slice_head(n = 25)
# 
# best
# 
# 
# 






