

## nikos AML data
library(DECIPHER)
library(Biostrings)
library(tidyverse)

source("/g/schwab/marco/repos/sabeRprobes/R/fncts.R")
otu_counts <- read_tsv("/g/schwab/marco/projects/osFISH/probe_design/nanopre_aml/2026_05_Nanopore_Share/2026_05_Nanopore_Share/final_output/otu_table.tsv")

otu_annotation <- read_tsv("/g/schwab/marco/projects/osFISH/probe_design/nanopre_aml/2026_05_Nanopore_Share/2026_05_Nanopore_Share/final_output/taxonomy_table.tsv")

otu_seq <- readDNAStringSet("/g/schwab/marco/projects/osFISH/probe_design/nanopre_aml/2026_05_Nanopore_Share/2026_05_Nanopore_Share/final_output/otus.fasta")


relevant <- c("Prorocentrum_micans", "Heterocapsa_rotundata", "Amphidinium_carterae")


barcode_ids_df <- tibble(
  species=relevant,
  barcode=c("FBF64117bc05","FBF64117bc04","FBF64117bc03")
)


rel_otus <- 
  otu_annotation %>%
  filter(Species %in% relevant) %>%
  pull(OTU)


ov_df <- otu_counts %>%
  left_join(otu_annotation %>% select(OTU, species=Species)) %>%
  filter(OTU %in% rel_otus) 



best_otus <-lapply(relevant, function(SPEC){
  
  print(SPEC)
  
  rel_bc <- barcode_ids_df %>% filter(species==SPEC) %>% pull(barcode)
  
  otu_best <- ov_df %>%
    filter(species==SPEC) %>%
    select(OTU,species, bc=all_of(rel_bc)) %>%
    filter(bc==max(bc)) %>%
    pull(OTU)
  
  best_seq <- otu_seq[which(names(otu_seq)==otu_best)]
  names(best_seq) <- SPEC#paste(SPEC, names(best_seq), sep="_")
  
  return(best_seq)
  
}) %>%
  Reduce(function(x,y){c(x,y)},.)


main_colors=c("green", "yellow","orange", "red","purple", "blue")

nOTU <- 1

all_probes <- lapply(seq(1, length(best_otus)), function(nOTU){
  
  
  
  rel_otu <- best_otus[nOTU]
  name <- names(rel_otu)
  
  print(name)
  
  all_possible_probes <- mine_probes_from_sequence(as.character(rel_otu), 41, 38, 3) %>%
    #rowwise() %>%
    mutate(
      len=nchar(seq),
      
      gc_content=str_count(seq, "G|C")/len
    ) %>%
    rowwise() %>%
    mutate(
      #tm=calc_tm_nn(seq),
      md_tm=get_oligo_tm(seq)
    )
  
  
  rel_temp <- all_possible_probes %>% 
    ## coarse filtering alredy here
    filter(
      between(gc_content, 0.35, 0.65),
      md_tm<75
    )%>%
    ungroup() %>%
    mutate(
      probe_id = row_number(),
      end = start + len - 1
    ) %>% 
    mutate(
      probe_organism=name
    )
  
  
  
  mapped <- rel_temp %>% select(probe_seq=seq, probe_id)  %>%
    crossing(target_id = names(best_otus)) %>%
    mutate(
      target_seq = as.character(best_otus[target_id]),
      
      aln = map2(
        probe_seq,
        target_seq,
        ~ pairwiseAlignment(
          pattern = DNAString(.x),
          subject = DNAString(.y),
          type = "local",
          substitutionMatrix = nucleotideSubstitutionMatrix(
            match = 1,
            mismatch = -1,
            baseOnly = TRUE
          ),
          gapOpening = -5,
          gapExtension = -2
        )
      ),
      
      score = map_dbl(aln, score),
      pid = map_dbl(aln, pid),
      nmatch = map_int(aln, nmatch),
      mismatch = map_int(aln, nmismatch),
      gaps = map_int(aln, ~ length(nindel(.x))),
      
      probe_start = map_int(aln, ~ start(pattern(.x))),
      probe_end   = map_int(aln, ~ end(pattern(.x))),
      
      target_start = map_int(aln, ~ start(subject(.x))),
      target_end   = map_int(aln, ~ end(subject(.x))),
      
      aligned_probe = map_chr(aln, ~ as.character(pattern(.x))),
      aligned_target = map_chr(aln, ~ as.character(subject(.x)))
    ) %>%
    select(
      probe_id,
      target_id,
      score,
      pid,
      nmatch,
      mismatch,
      gaps,
      probe_start,
      probe_end,
      target_start,
      target_end,
      aligned_probe,
      aligned_target
    ) %>%
    left_join(rel_temp, by="probe_id")
  
  
  unspecifics <- mapped %>%
    filter(target_id!=probe_organism) %>% 
    mutate(unaligned_len=1-nmatch/len) %>% 
    group_by(probe_id) %>%
    summarize(
      product=prod(unaligned_len),
      summed=sum(unaligned_len),
      minimal=min(unaligned_len)
    )
  
  
  
  joined <- unspecifics %>%
    filter(minimal>0.1) %>%
    left_join(rel_temp, by="probe_id") 
  
  
  grouped <- lapply(unique(joined$md_tm), function(TMP){
    
    joined %>%
      filter(md_tm==TMP) %>%
      assign_nonoverlap_groups_by_product(., 4) %>%
      mutate(
        temp_group=TMP
        #group_id=paste(TMP, nonoverlap_group, sep="_")
      )
    
  }) %>%
    bind_rows() %>%
    left_join(
      group_by(., temp_group, nonoverlap_group) %>% tally()
    ) %>%
    #filter(n>2) %>%
    mutate_at(.vars="nonoverlap_group", .funs=as.character)
  
  return(grouped) 
  
})


# 
# ggplot(grouped %>%
#          filter(n>2), aes(x=start, xend=end, y=nonoverlap_group, yend=nonoverlap_group, color=product))+
#   geom_segment(size=5)+
#   facet_grid(temp_group+probe_organism~"1", scales="free_y", space="free")+
#   scale_color_gradientn(colors=main_colors)
# 


pd <- bind_rows(all_probes) %>%
  filter(n>2) %>%
  group_by(probe_organism, temp_group, nonoverlap_group) %>%
  arrange(desc(product)) %>%
  mutate(int_group_rank=c(1:length(product))) %>%
  #filter(int_group_rank<=3) %>%
  left_join(
    .,
    summarize(., mn_prod=mean(product)) %>% arrange(desc(mn_prod))%>%
   group_by(probe_organism, temp_group)%>%
   mutate(top5_rank=c(1:length(mn_prod)))
  )  %>%
   arrange(probe_organism, temp_group)  %>%
   filter(top5_rank<=5)
  


final_plot <- ggplot(pd, aes(x=start, xend=end, y=as.character(top5_rank), yend=as.character(top5_rank), color=product))+
  geom_segment(size=5)+
  geom_text(aes(label=probe_id, x=start+18), color="black")+
  facet_grid(temp_group~probe_organism, scales="free_y", space="free")+
  scale_color_gradientn(colors=main_colors)+
  theme_bw()





outdir="/home/rheinnec"
plot_name <- "test_full"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=20, height=6)
final_plot
dev.off()


## final probe selections


selected <- tibble::tribble(
  ~organism, ~probe_id,   
  "Amphidinium_carterae", 86,
  "Amphidinium_carterae", 228,
  "Amphidinium_carterae", 310,
  "Heterocapsa_rotundata", 18,
  "Heterocapsa_rotundata", 20,
  "Heterocapsa_rotundata", 32,
  "Prorocentrum_micans", 30,
  "Prorocentrum_micans", 54,
  "Prorocentrum_micans", 76,
  
)



final_probes <- bind_rows(all_probes) %>%
  mutate(full_id=paste(probe_organism, probe_id, sep="_")) %>%
  filter(
    md_tm==74,
    full_id %in% paste(selected$organism, selected$probe_id, sep="_")
  )



msa <- DECIPHER::AlignSeqs(best_otus)


######### from HERE CHATGP test

msa_df <- tibble(
  probe_organism = names(msa),
  aln_seq = as.character(msa)
) %>%
  mutate(bases = strsplit(aln_seq, "")) %>%
  select(probe_organism, bases) %>%
  unnest_longer(bases, indices_to = "aln_pos") %>%
  group_by(probe_organism) %>%
  mutate(
    ungapped_pos = cumsum(bases != "-"),
    ungapped_pos = if_else(bases == "-", NA_integer_, ungapped_pos)
  ) %>%
  ungroup()

probe_ranges <- final_probes %>%
  select(
    full_id,
    probe_organism,
    start,
    end,
    nonoverlap_group,
    temp_group,
    product
  )


probe_ranges_msa <- final_probes %>%
  select(full_id, probe_organism, start, end) %>%
  left_join(
    msa_df %>%
      filter(!is.na(ungapped_pos)) %>%
      select(probe_organism, ungapped_pos, aln_pos),
    by = "probe_organism",
    relationship = "many-to-many"
  ) %>%
  filter(ungapped_pos >= start, ungapped_pos <= end) %>%
  group_by(full_id, probe_organism) %>%
  summarise(
    probe_aln_start = min(aln_pos),
    probe_aln_end = max(aln_pos),
    .groups = "drop"
  )
msa_plot_df2 <- msa_df %>%
  left_join(probe_ranges_msa, by = "probe_organism", relationship = "many-to-many") %>%
  mutate(
    in_own_probe = aln_pos >= probe_aln_start & aln_pos <= probe_aln_end
  ) %>%
  group_by(probe_organism, aln_pos, bases, ungapped_pos) %>%
  summarise(
    in_own_probe = any(in_own_probe, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(aln_pos) %>%
  mutate(
    column_has_probe = any(in_own_probe),
    probe_base = bases[in_own_probe & bases != "-"][1],
    
    is_probe_row = in_own_probe,
    
    is_nonprobe_mismatch =
      column_has_probe &
      !is_probe_row &
      bases != "-" &
      !is.na(probe_base) &
      bases != probe_base,
    
    is_nonprobe_gap =
      column_has_probe &
      !is_probe_row &
      bases == "-",
    
    fill_status = case_when(
      is_probe_row & bases == "A" ~ "probe_A",
      is_probe_row & bases == "C" ~ "probe_C",
      is_probe_row & bases == "G" ~ "probe_G",
      is_probe_row & bases == "T" ~ "probe_T",
      is_probe_row & bases == "-" ~ "gap_in_probe",
      
      is_nonprobe_mismatch & bases == "A" ~ "mismatch_A",
      is_nonprobe_mismatch & bases == "C" ~ "mismatch_C",
      is_nonprobe_mismatch & bases == "G" ~ "mismatch_G",
      is_nonprobe_mismatch & bases == "T" ~ "mismatch_T",
      
      is_nonprobe_gap ~ "gap_in_other_seq",
      
      TRUE ~ "background"
    )
  ) %>%
  ungroup()

padding <- 20
merge_distance <- 30  # merge windows if closer than this

probe_windows <- msa_plot_df2 %>%
  filter(in_own_probe) %>%
  distinct(aln_pos) %>%
  arrange(aln_pos) %>%
  mutate(
    new_window = aln_pos > lag(aln_pos, default = first(aln_pos)) + merge_distance,
    window_id = cumsum(new_window) + 1
  ) %>%
  group_by(window_id) %>%
  summarise(
    window_start = max(1, min(aln_pos) - padding),
    window_end   = max(aln_pos) + padding,
    .groups = "drop"
  )
msa_plot_df_zoom <- msa_plot_df2 %>%
  crossing(probe_windows) %>%
  filter(aln_pos >= window_start, aln_pos <= window_end) %>%
  mutate(
    window_label = paste0("MSA ", window_start, "-", window_end),
    fill_status=case_when(
      bases=="-"&is.na(probe_base) ~ "gap_in_probe",
      is_probe_row&!is_nonprobe_mismatch ~ "probe_aligned",
      TRUE ~ fill_status
    )
      
  )

tp <- ggplot(msa_plot_df_zoom, aes(x = aln_pos, y = probe_organism)) +
  geom_tile(aes(fill = fill_status, alpha=is_probe_row), color = "white", linewidth = 0.08) +
  geom_text(aes(label = bases), size = 2.2, family = "mono") +
  facet_wrap(~ window_label, scales = "free_x", ncol = 1) +
  scale_alpha_manual(values=c(0.5, 1))+
  scale_fill_manual(
    values = c(
      "background" = "grey90",
      "probe_aligned" = "grey70",
      "probe_A" = "#4DAF4A",
      "probe_C" = "#377EB8",
      "probe_G" = "#E41A1C",
      "probe_T" = "#984EA3",
      
      "mismatch_A" = "#4DAF4A",
      "mismatch_C" = "#377EB8",
      "mismatch_G" = "#E41A1C",
      "mismatch_T" = "#984EA3",
      
      "gap_in_probe" = "#E66100",
      "gap_in_other_seq" = "#E66100"
    ),
    name = "Position type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 9),
    strip.background = element_rect(fill = "grey95")
  ) +
  labs(
    x = "MSA position",
    y = NULL
  )


outdir="/home/rheinnec"
plot_name <- "test_final"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=10, height=7)
tp
dev.off()




primers <- tibble(
  
  primer_id=c(25, 27, 28),
  #pool_id=c("p1", "p2", "p3"),
  primer_seq=c("CCAATAATA",
               "CATCATCAT",
               "CAACTTAAC"),
  id_primer=c("p25", "p27", "p28")
  
)


dna <- DNAStringSet(
  setNames(final_probes$seq, final_probes$full_id)
) %>% reverseComplement() 


finalized_order <-  
  tibble(
    full_id = names(dna),
    probe_seq = as.character(dna)
  ) %>%
  mutate(
    name=str_replace(full_id, "Prorocentrum_micans", "PMp28") %>%
      str_replace("Heterocapsa_rotundata", "HRp27") %>%
      str_replace("Amphidinium_carterae", "ACp25"),
    id_primer=str_extract(name, "p25|p27|p28")
  ) %>%
  left_join(primers, by="id_primer")%>%
  mutate(
    seq_full=
      paste0(
        ## target sequence
        probe_seq,
        ## linker sequence
        "TTTTT",
        "A",
        primer_seq,
        "A",
        primer_seq
      )
  ) %>%
    rowwise() %>%
    mutate(
      order_name=str_remove(name, id_primer) %>% paste(.,id_primer, sep="_")
    ) 
  
  
to_order <- finalized_order %>%
  select(order_name, seq_full)
  
write_tsv(finalized_order, file="/g/schwab/marco/projects/osFISH/orders/order_26-05-12.tsv")    



  
  final_probe_sequences <- to_order %>%
  ## adding scrambled sequence... no hits in BLAST and PR2 database
  # bind_rows(tibble(id="scrambled_N1", revCom="CTACGGTGGGACGAGAAATCTTACACAATCTGTGCGAAGT", 
  #                  id_primer=c("p27", "p28"))) %>%
  select(id, revCom, id_primer)%>%
  left_join(primers, by="id_primer") %>%
  rowwise() ,  
  





##############

###################
#########################
######################################



pair_df <- expand_grid(
  otu1 = names(best_otus),
  otu2 = names(best_otus)
) %>%
  filter(otu1 < otu2)


alignments <- lapply(seq(1, nrow(pair_df)), function(nP){
  
  
  seq1 <- best_otus[pair_df[[nP, "otu1"]]]
  seq2 <- best_otus[pair_df[[nP, "otu2"]]]
  
  aln <- pairwiseAlignment(seq1, seq2, type = "global")
  
  return(aln)
  
})




diffs <- pairwise_difference_vectors(alignments)




diffs_long <- map_dfr(diffs, pairwise_to_long)



per_sequence_position <- diffs_long %>%
  group_by(seq_name, seq_pos) %>%
  summarise(
    base = dplyr::first(base),
    n_comparisons = n(),
    n_different = sum(different),
    n_gaps_vs_others = sum(gap),
    n_mismatches_vs_others = sum(mismatch),
    total_score = sum(score),
    different_to = paste(other_seq[different], collapse = "; "),
    .groups = "drop"
  )




windows <- per_sequence_position %>%
  group_by(seq_name) %>%
  group_split() %>%
  map_dfr(sliding_windows_per_seq, window = 40) %>%
  ungroup() %>%
  mutate(tot_norm=4*total_score/(max(.$total_score)))

p <- ggplot(per_sequence_position, aes(x=seq_pos, y=total_score))+
  facet_wrap(~seq_name, ncol=1)+
  geom_line()+
  geom_text(aes(label=base, y=-1), size=1.5)+
  geom_line(data=windows, aes(x=start_pos, y=tot_norm), color="red")+
  geom_col(data=windows, aes(x=start_pos, y=tot_norm, fill=total_score), alpha=0.5)+
  scale_fill_gradientn(colors = main_colors)+
  scale_x_continuous(breaks=seq(1, max(per_sequence_position$seq_pos), 10))+
  theme_bw()



outdir="/home/rheinnec"
plot_name <- "test"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=70, height=5)
p
dev.off()




mapped    %>% View()





group_sequences(rel_temp, 5)


probes <- 
  tibble(
    ac1 = "AGTCATAGAGAATCGGCGCAGGCTCTGCCTTGCTA",
    ac2 = "TTACAACAGCAATAATTTTCGCAGTGCTTCATCGCACATGGAT",
    ac3 = "TGGCTAAGTAGCTATGCATAGCTGTGGTTTTGTGGACAGT",
    
    hr1 = "AGTGTTCGGCAACGAGCGTTGCAGCGGAAAGTTTAG",  #"CAGCAGTGTTCGGCAACGAGCGTTGCAGCGGAAAGTTTAG",
    hr2 = "GAATGAGTAGAATTTAAAACCCTTTACGAGTACCGATTGGAGGG", 
    hr3 = "TATTTGATGGTCATTCTTACATGGATAACCGTGGTAATTCTAGAGC",   # "TAATTCGGACTGCAGCAGTGTTCGGCAACGAGCGTTG",
    hr4 = "GCCAAAACCCGACTTCTTGGAAGGGTTGTGTTTATTAGTTAC",
    
    pm1 = "TGATGGTCACTCTTTACATGGATAACTGTGCTAATTGTAGAGC" , #"TTATTTGATGGTCACTCTTTACATGGATAACTGTGCTAATTGTAGAGCTAA" 
    pm2 = "CGGAGAGGGAGAAACTCTGTCTGAGAAACGGCTACCACAT",
    pm3 = "TAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCA"
  ) %>%
  gather() %>%
  rowwise() %>%
  mutate(
    len=nchar(value),
    tm=calc_tm_nn(value),
    gc_content=str_count(value, "G|C")/len
  ) %>%
  arrange(tm)


probes

patt <- "etero"
seq <- probes %>%
  filter(key=="hr3") %>%
  pull(value)


pairwiseAlignment(seq, best_otus[which(str_detect(names(best_otus), patt))], type = "local")




probe_align_scores <- probes %>%
  select(probe_id = key, probe_seq = value) %>%
  crossing(target_id = names(best_otus)) %>%
  mutate(
    target_seq = as.character(best_otus[target_id]),
    
    aln = map2(
      probe_seq,
      target_seq,
      ~ pairwiseAlignment(
        pattern = DNAString(.x),
        subject = DNAString(.y),
        type = "local",
        substitutionMatrix = nucleotideSubstitutionMatrix(
          match = 1,
          mismatch = -1,
          baseOnly = TRUE
        ),
        gapOpening = -5,
        gapExtension = -2
      )
    ),
    
    score = map_dbl(aln, score),
    pid = map_dbl(aln, pid),
    nmatch = map_int(aln, nmatch),
    mismatch = map_int(aln, nmismatch),
    gaps = map_int(aln, ~ length(nindel(.x))),
    
    probe_start = map_int(aln, ~ start(pattern(.x))),
    probe_end   = map_int(aln, ~ end(pattern(.x))),
    
    target_start = map_int(aln, ~ start(subject(.x))),
    target_end   = map_int(aln, ~ end(subject(.x))),
    
    aligned_probe = map_chr(aln, ~ as.character(pattern(.x))),
    aligned_target = map_chr(aln, ~ as.character(subject(.x)))
  ) %>%
  select(
    probe_id,
    target_id,
    score,
    pid,
    nmatch,
    mismatch,
    gaps,
    probe_start,
    probe_end,
    target_start,
    target_end,
    aligned_probe,
    aligned_target
  )



probe_align_scores %>% left_join(probes, by=c("probe_id"="key")) %>% mutate(aligned_len=nmatch/len) %>% arrange(desc(aligned_len)) %>% View()




lapply(seq(1, nrow(pot_probes)), function(nP){
  
  
  
  
})


#calc_tm_nn("TTACCCGTCATTGCCACGGTAAGCCCATATCCTACCTTTATAATA")


patt <- "etero"
seq <- "TTTTATACGACGAAACTGCGAATGGCTCATTAAAACAGT"


pairwiseAlignment(seq, best_otus[which(str_detect(names(best_otus), patt))], type = "local")

n <- 68

best_otus[which(str_detect(names(best_otus), patt))] %>% as.character() %>% str_sub(n, n+50)



difference_vectors <- per_sequence_position %>%
  arrange(seq_name, seq_pos) %>%
  group_by(seq_name) %>%
  summarise(
    diff_vector = list(total_score),
    .groups = "drop"
  )












best_seq <- otu_seq[which(names(otu_seq) %in% best_otus$OTU)]



hr_probe <- reverseComplement(DNAString("AATACCGCACCACACAGTCAAGTGCAGATACGTTCTCCAA"))
hr_seq <- otu_seq["cluster2266_FBF64117bc04"]


alignment <- pairwiseAlignment(hr_probe, hr_seq, type = "local")
print(alignment)


pm_probe <- reverseComplement(DNAString("CTCGAAGTCGGGTTTGGGCGCATGTATTAGCTCTACAAT"))
pm_seq <- otu_seq["cluster1383_FBF64117bc05"]


alignment <- pairwiseAlignment(hr_probe, pm_seq, type = "local")
print(alignment)

writePairwiseAlignments(alignment)



ac_probe=reverseComplement(DNAString("TTTGCTGCACCCTTCCTCAGCACTTGAGCCTAGAGAT"))

ac_seq <- otu_seq["cluster16365_FBF64117bc07"]


alignment <- pairwiseAlignment(ac_probe, ac_seq, type = "local")
print(alignment)

writePairwiseAlignments(alignment)


score_msa_variability <- function(msa, window = 20, step = 1) {
  
  # convert DECIPHER alignment to character matrix
  seqs <- as.character(msa)
  aln_mat <- do.call(rbind, strsplit(seqs, split = ""))
  
  n_seq <- nrow(aln_mat)
  aln_len <- ncol(aln_mat)
  
  # per-column score
  col_scores <- lapply(seq_len(aln_len), function(i) {
    x <- aln_mat[, i]
    
    n_gaps <- sum(x == "-")
    
    # pairwise comparisons
    pairs <- combn(seq_len(n_seq), 2)
    pair_vals <- apply(pairs, 2, function(p) {
      a <- x[p[1]]
      b <- x[p[2]]
      
      mismatch <- a != b
      gap_in_pair <- a == "-" || b == "-"
      
      data.frame(
        mismatch = mismatch,
        gap_pair = gap_in_pair
      )
    })
    
    pair_vals <- do.call(rbind, pair_vals)
    
    data.frame(
      aln_pos = i,
      n_gaps = n_gaps,
      n_pairwise_mismatches = sum(pair_vals$mismatch),
      n_gap_pairs = sum(pair_vals$gap_pair),
      variability_score = sum(pair_vals$mismatch) + n_gaps
    )
  })
  
  col_scores <- do.call(rbind, col_scores)
  
  # sliding-window summary
  starts <- seq(1, aln_len - window + 1, by = step)
  
  windows <- lapply(starts, function(s) {
    e <- s + window - 1
    sub <- col_scores[s:e, ]
    
    data.frame(
      start = s,
      end = e,
      window_size = window,
      total_gaps = sum(sub$n_gaps),
      total_pairwise_mismatches = sum(sub$n_pairwise_mismatches),
      total_gap_pairs = sum(sub$n_gap_pairs),
      variability_score = sum(sub$variability_score)
    )
  })
  
  windows <- do.call(rbind, windows)
  windows <- windows[order(-windows$variability_score), ]
  
  list(
    per_column = col_scores,
    windows = windows
  )
}



pre_msa_list1 <- c(pm_seq, hr_seq)
pre_msa_list2 <- c(pm_seq, ac_seq)
pre_msa_list3 <- c(ac_seq, hr_seq)


tot <- lapply(list(pre_msa_list1, pre_msa_list2, pre_msa_list3), function(ALN){
  
  
  msa <- DECIPHER::AlignSeqs(ALN)
  
  res <- score_msa_variability(msa, window = 30, step = 1) 
  
  return(res$per_column)
  
})


dat <- lapply(1:length(tot), function(N){
  tot[[N]] %>% mutate(comp=N)
}) %>%
  bind_rows() %>%
  as_tibble()



p <- ggplot(dat, aes(x=aln_pos, y=variability_score))+
  facet_wrap(~comp, ncol=1)+
  geom_line()+
  scale_x_continuous(breaks=seq(1, max(dat$aln_pos), 10))



outdir="/home/rheinnec"
plot_name <- "test"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=70, height=5)
p
dev.off()






head(res$windows, 10)

DECIPHER::BrowseSeqs(msa)


#











#1: TRECC 16 C. Neogracilis
#2: TRECC 108 bentic diatom
#3: TRECC 134 C. peruvienus
#4: Coscinodiscus granii
#5: Pyrecystis iunula
#6: Odontella sinensis
#7: RCC 1212 enut
#8: AC Amphidinium carterae (no thecal plates)
#9: HR Heterocapsa rotundata (thecal plates, small ~ 10 um)
#10: PM Prorocentrum micans (thick thecal plates, big ~ 40 um)




#2 → (BC-01)
#6 → (BC-02)
#8 → (BC-03)
#9 → (BC-04)
#10 → (BC-05)
#mix of the 5 above equimol → (BC-06)
#mix of the 5 above equivolume → (BC-07)






library(ShortRead)


















fq <- readFastq("/scratch/rheinnec/taxseq/re_basecalling/fastqs/fastq_pass/barcode04/FBF64117_pass_barcode04_16e02936_fc9fa64b_0.fastq")

writeFastq(fq[1:1000], file = "/scratch/rheinnec/taxseq/filtered_barco04_1000.fastq")

write_tsv(ref2taxid_file, file="/scratch/rheinnec/taxseq/pr2_ref2taxid.tsv")

# 
# fq
# 
sread(fq)[1:20]


msa <- DECIPHER::AlignSeqs(sread(fq)[1:20])

DECIPHER::BrowseSeqs(sread(fq)[1:20])



saveRDS(msa, file="/scratch/rheinnec/taxseq/barcode04_msa_full.rds")

library(Biostrings)

opt <- tibble(
  #  file_target="/g/schwab/marco/projects/osFISH/test2/Prorocentrum_micans_seq_all.fasta",
  file_reference="/scratch/rheinnec/np_test/pr2_version_5.1.1_SSU_taxo_long.fasta",
  # temp_min=55,
  # temp_max=75,
  # length=40,
  # length_range=3
)


ref_seq <- readDNAStringSet(opt$file_reference)

ref_seq_corr_names <- ref_seq

names(ref_seq_corr_names) <- str_remove_all(names(ref_seq_corr_names), "\"")



ref2taxid_file <- tibble(
  refname=names(ref_seq_corr_names)
) %>%
  mutate(taxid=c(1:nrow(.)))

write_tsv(ref2taxid_file, file="/scratch/rheinnec/np_test/reftax.tsv", col_names = F)


writeXStringSet(ref_seq_corr_names, file="/scratch/rheinnec/np_test/pr2_namecorr.fasta")






HR <- ref_seq[which(str_detect(names(ref_seq), "eterocapsa"))]

hr_long <- HR[which(width(HR)>3000)]


barcode04 <- sread(fq)[1:20]

c(hr_long, barcode04)



seq1 <-  hr_long[[1]]
seq2 <-  barcode04[[1]]

alignment <- pairwiseAlignment(seq1, seq2, type = "global")
DECIPHER::BrowseSeqs(alignment)


msa <- DECIPHER::AlignSeqs(c(hr_long, barcode04[1]))








