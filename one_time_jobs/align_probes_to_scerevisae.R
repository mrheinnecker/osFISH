
library(googlesheets4)
library(tidyverse)
library(Biostrings)

## function def


# x = list of LocalPairwiseAlignmentsSingleSubject objects
extract_hits <- function(x) {
  do.call(rbind, lapply(seq_along(x), function(i) {
    a <- x[[i]]
    
    # aligned strings (include gaps)
    pat_aln <- as.character(pattern(a))
    sub_aln <- as.character(subject(a))
    
    # ungapped lengths (true matched span lengths in each sequence)
    pat_ungapped <- nchar(gsub("-", "", pat_aln))
    sub_ungapped <- nchar(gsub("-", "", sub_aln))
    
    # start positions reported by Biostrings (1-based)
    pat_start <- start(pattern(a))
    sub_start <- start(subject(a))
    
    # end positions (1-based, inclusive)
    pat_end <- pat_start + pat_ungapped - 1
    sub_end <- sub_start + sub_ungapped - 1
    
    data.frame(
      hit_id      = i,
      pattern_start= pat_start,
      pattern_end  = pat_end,
      subject_start= sub_start,
      subject_end  = sub_end,
      pattern_ungapped = pat_ungapped,
      subject_ungapped  = sub_ungapped,
      score       = score(a),
      pattern_aln = pat_aln,
      subject_aln = sub_aln,
      stringsAsFactors = FALSE
    )
  }))
}



# aln_list: named list of LocalPairwiseAlignmentsSingleSubject
# target_sequences: tibble with probe_id, target_seq
# (these are from your earlier code)

alignment_map <- function(a, hit_id) {
  pat <- strsplit(as.character(pattern(a)), "", fixed = TRUE)[[1]]
  sub <- strsplit(as.character(subject(a)), "", fixed = TRUE)[[1]]
  stopifnot(length(pat) == length(sub))
  
  s_pos <- start(subject(a)) - 1
  p_pos <- start(pattern(a)) - 1
  
  out <- vector("list", length(pat))
  for (k in seq_along(pat)) {
    if (sub[k] != "-") s_pos <- s_pos + 1
    if (pat[k] != "-") p_pos <- p_pos + 1
    
    out[[k]] <- tibble(
      hit_id = hit_id,
      aln_col = k,
      pattern_pos = if (pat[k] == "-") NA_integer_ else p_pos,
      subject_pos = if (sub[k] == "-") NA_integer_ else s_pos,
      pat_base = pat[k],
      sub_base = sub[k],
      is_gap = (pat[k] == "-" | sub[k] == "-"),
      is_mismatch = (! (pat[k] == "-" | sub[k] == "-") & pat[k] != sub[k])
    )
  }
  bind_rows(out)
}

# Make a tidy "letters plotted at subject coordinates" table for all probes
make_subject_letter_df <- function(aln_list, target_sequences_tbl) {
  map_dfr(seq_along(aln_list), function(i) {
    a <- aln_list[[i]]
    probe_id <- target_sequences_tbl$probe_id[i]
    full_target <- target_sequences_tbl$target_seq[i]
    
    # mapping within the aligned window (handles gaps)
    m <- alignment_map(a, hit_id = probe_id)
    
    # pull the non-gap pattern<->subject mapping
    core <- m %>%
      filter(!is.na(pattern_pos)) %>%
      group_by(pattern_pos) %>%
      summarise(
        subject_pos = subject_pos[which(!is.na(subject_pos))[1]],  # may be NA if aligned to subject gap
        is_mismatch = any(is_mismatch, na.rm = TRUE),
        .groups = "drop"
      )
    
    # aligned window bounds on the PROBE (pattern)
    pat_start <- start(pattern(a))
    pat_end <- pat_start + nchar(gsub("-", "", as.character(pattern(a)))) - 1
    
    # subject window bounds
    sub_start <- start(subject(a))
    sub_end <- sub_start + nchar(gsub("-", "", as.character(subject(a)))) - 1
    
    # build full probe base table
    full_df <- tibble(
      probe_id = probe_id,
      probe_pos = seq_len(nchar(full_target)),
      base = strsplit(full_target, "", fixed = TRUE)[[1]]
    ) %>%
      left_join(core, by = c("probe_pos" = "pattern_pos")) %>%
      mutate(
        # fill in subject_pos for unaligned flanks by extrapolation
        subject_pos = case_when(
          !is.na(subject_pos) ~ subject_pos,
          
          # left flank (before aligned start on probe)
          probe_pos < pat_start ~ sub_start - (pat_start - probe_pos),
          
          # right flank (after aligned end on probe)
          probe_pos > pat_end ~ sub_end + (probe_pos - pat_end),
          
          # inside aligned window but aligned to subject gap => keep NA (won't plot)
          TRUE ~ NA_integer_
        ),
        is_aligned = probe_pos >= pat_start & probe_pos <= pat_end
      ) %>%
      # attach useful bounds for zooming/labeling
      mutate(
        subject_start = sub_start,
        subject_end = sub_end
      )
    
    full_df
  })
}

pad_to_full_target <- function(full_target, pat_start, pat_end, aligned_pattern_gapped) {
  # aligned_pattern_gapped is the gapped aligned window (e.g., A-C-T)
  left_pad  <- if (pat_start > 1) paste(rep(" ", pat_start - 1), collapse = "") else ""
  right_pad <- if (pat_end < nchar(full_target)) paste(rep(" ", nchar(full_target) - pat_end), collapse = "") else ""
  paste0(left_pad, aligned_pattern_gapped, right_pad)
}



extract_hits_with_full_target <- function(aln_list, target_sequences_tbl) {
  map_dfr(seq_along(aln_list), function(i) {
    a <- aln_list[[i]]
    probe_id <- target_sequences_tbl$probe_id[i]
    full_target <- target_sequences_tbl$target_seq[i]
    
    pat_aln <- as.character(pattern(a))   # aligned window (gapped)
    sub_aln <- as.character(subject(a))   # aligned window (gapped)
    
    # start on full target (pattern)
    pat_start <- start(pattern(a))
    
    # aligned span on the target WITHOUT gaps
    pat_ungapped <- nchar(gsub("-", "", pat_aln))
    pat_end <- pat_start + pat_ungapped - 1
    
    # start/end on the long subject
    sub_start <- start(subject(a))
    sub_ungapped <- nchar(gsub("-", "", sub_aln))
    sub_end <- sub_start + sub_ungapped - 1
    
    # make the *full* target visible, highlighting the aligned region
    before <- if (pat_start > 1) substr(full_target, 1, pat_start - 1) else ""
    mid    <- substr(full_target, pat_start, pat_end)
    after  <- if (pat_end < nchar(full_target)) substr(full_target, pat_end + 1, nchar(full_target)) else ""
    
    target_seq_marked <- paste0(before, "[", mid, "]", after)
    
    tibble(
      probe_id = probe_id,
      score = score(a),
      pattern_start = pat_start,
      pattern_end = pat_end,
      subject_start = sub_start,
      subject_end = sub_end,
      aligned_pattern = pat_aln,
      aligned_subject = sub_aln,
      target_seq_full = full_target,
      target_seq_marked = target_seq_marked
    )
  })
}







## functiond ef done

probe_ov_url <- "https://docs.google.com/spreadsheets/d/1vNJssytzfJEsmYrIdLDZi3v_V6m-sokpKpTKtSpLXbk/edit?gid=1469504260#gid=1469504260"

probe_ov <- probe_ov_url %>% read_sheet(sheet="probes_info", col_types="c") 

soi <- c("AS2p28","PL2p28", "TH3p27","AC2p27","HR3p27","P4p27")

target_sequences <- probe_ov %>%
  select(probe_id, seq=Sequence) %>%
  filter(probe_id %in% soi) %>%
  mutate(target_seq_raw=str_remove(seq, "TTTTTACATCATCATACATCATCAT$|TTTTTACAACTTAACACAACTTAAC")) %>%
  mutate(target_seq = target_seq_raw %>% 
           DNAStringSet() %>% reverseComplement() %>% as.character()) %>%
  select(-seq, -target_seq_raw)


rRNA_s_cerevisiae="TATCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCAATTTATACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATAGTTCCTTTACTACATGGTATAACTGTGGTAATTCTAGAGCTAATACATGCTTAAAATCTCGACCCTTTGGAAGAGATGTATTTATTAGATAAAAAATCAATGTCTTCGGACTCTTTGATGATTCATAATAACTTTTCGAATCGCATGGCCTTGTGCTGGCGATGGTTCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTTTCAACGGGTAACGGGGAATAAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTAATTCAGGGAGGTAGTGACAATAAATAACGATACAGGGCCCATTCGGGTCTTGTAATTGGAATGAGTACAATGTAAATACCTTAACGAGGAACAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACTTTGGGCCCGGTTGGCCGGTCCGATTTTTTCGTGTACTGGATTTCCAACGGGGCCTTTCCTTCTGGCTAACCTTGAGTCCTTGTGGCTCTTGGCGAACCAGGACTTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCGTATTGCTCGAATATATTAGCATGGAATAATAGAATAGGACGTTTGGTTCTATTTTGTTGGTTTCTAGGACCATCGTAATGATTAATAGGGACGGTCGGGGGCATCAGTATTCAATTGTCAGAGGTGAAATTCTTGGATTTATTGAAGACTAACTACTGCGAAAGCATTTGCCAAGGACGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCGAAGATGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCGGGTGGTGTTTTTTTAATGACCCACTCGGCACCTTACGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTCACCAGGTCCAGACACAATAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTTGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGCTTAATTGCGATAACGAACGAGACCTTAACCTACTAAATAGTGGTGCTAGCATTTGCTGGTTATCCACTTCTTAGAGGGACTATCGGTTTCAAGCCGATGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGACGTTCTGGGCCGCACGCGCGCTACACTGACGGAGCCAGCGAGTCTAACCTTGGCCGAGAGGTCTTGGTAATCTTGTGAAACTCCGTCGTGCTGGGGATAGAGCATTGTAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTAGTACCGATTGAATGGCTTAGTGAGGCCTCAGGATCTGCTTAGAGAAGGGGGCAACTCCATCTCAGAGCGGAGAATTTGGACAAACTTGGTCATTTAGAGGAACTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA"


subject <- DNAString(rRNA_s_cerevisiae)

sm <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)



aln_list <- pmap(
  target_sequences,
  function(probe_id, target_seq) {
    pairwiseAlignment(
      DNAString(target_seq), 
      subject, 
      type = "local",
      substitutionMatrix = sm,
      gapOpening = -4,      # less negative => easier to open a gap
      gapExtension = -1     # less negative => easier to extend a gap
      
      )
  }
)
names(aln_list) <- target_sequences$probe_id




hits_df <- hits_df %>%
  mutate(
    pattern_aln_padded = pmap_chr(
      list(target_seq_full, pattern_start, pattern_end, aligned_pattern),
      pad_to_full_target
    )
  )

hits_df$pattern_aln_padded[1]






subject_len <- length(subject)

hits_df <- hits_df %>%
  mutate(probe_id = factor(probe_id, levels = rev(probe_id)))


# hits_df must contain: probe_id, target_seq_full, pattern_start, pattern_end
# (from the previous extract_hits_with_full_target() step)

letters_df <- hits_df %>%
  transmute(
    probe_id,
    target_seq = target_seq_full,
    pattern_start,
    pattern_end
  ) %>%
  mutate(
    bases = str_split(target_seq, "", simplify = FALSE),
    pos   = map(bases, seq_along)
  ) %>%
  unnest(c(bases, pos)) %>%
  dplyr::rename(base = bases) %>%
  mutate(
    is_aligned = pos >= pattern_start & pos <= pattern_end,
    probe_id = factor(probe_id, levels = rev(unique(probe_id)))
  )



letter_df <- make_subject_letter_df(aln_list, target_sequences) %>%
  mutate(probe_id = factor(probe_id, levels = rev(unique(probe_id))))

# Optional: choose a zoom window around all hits (recommended; otherwise it can get very wide)
xmin <- min(letter_df$subject_pos, na.rm = TRUE) - 10
xmax <- max(letter_df$subject_pos, na.rm = TRUE) + 10


scer_df <- tibble(
  subject_pos=seq(1,nchar(rRNA_s_cerevisiae)),
  base=  lapply(seq(1,nchar(rRNA_s_cerevisiae)), function(P){
    return(str_sub(rRNA_s_cerevisiae, P,P))
  }) %>% unlist()
)


lab_df <- letter_df %>%
  group_by(probe_id) %>%
  filter(
    probe_pos==max(probe_pos)
  ) %>%
  mutate(subject_pos=subject_pos+2)

p <- ggplot(letter_df %>% filter(!is.na(subject_pos)),
       aes(x = subject_pos, y = probe_id, label = base)) +

  geom_tile(aes(fill=base))+
  geom_text(aes(fontface = ifelse(is_aligned, "bold", "plain")),
            size = 4) + 
  geom_tile(data=scer_df,aes(fill=base, y="s.cerevisae"))+
  geom_text(data=scer_df, aes(y="s.cerevisae"))+
  geom_text(data=lab_df, aes(x=subject_pos, y=probe_id, label=probe_id), hjust=0)+
  
 
  coord_cartesian(xlim = c(xmin, xmax)) +
  labs(
    x = "rRNA_s_cerevisiae coordinate (subject)",
    y = NULL,
    title = "Full probe sequences placed on the subject coordinate axis",
    subtitle = "Aligned region in bold; flanks extrapolated"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,1800, 25))+
  theme(
    text = element_text(family = "mono"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


# 
# 
# 
# 
# 
# p <- ggplot(letter_df %>% filter(!is.na(subject_pos)),
#        aes(x = subject_pos, y = "t1", label = base)) +
#   geom_text(aes(fontface = ifelse(is_aligned, "bold", "plain"), color=probe_id),
#             size = 4)+ 
#   geom_point(data = letter_df %>% filter(!is.na(subject_pos) & is_mismatch),
#                                   aes(x = subject_pos, y = probe_id), size = 2) +
#   
#   ## rRNA scerevisae seqeunce
#   geom_text(data=scer_df, aes(y="t2"))+
#   
#   coord_cartesian(xlim = c(xmin, xmax)) +
#   labs(
#     x = "rRNA_s_cerevisiae coordinate (subject)",
#     y = NULL,
#     title = "Full probe sequences placed on the subject coordinate axis",
#     subtitle = "Aligned region in bold; flanks extrapolated"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(family = "mono"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )

pdf("/g/schwab/marco/seqov_gap4_ext1_mm3.pdf", width=200, height=2)
p
dev.off()


# 
# 
# all_alignments <- apply(target_sequences, 1, function(ROW){
#   
#   alignment <- pairwiseAlignment(DNAString(ROW[["target_seq"]]), 
#                                  DNAString(rRNA_s_cerevisiae), type = "local")
#   
#   
#   return(alignment)
#   
# })
# 


