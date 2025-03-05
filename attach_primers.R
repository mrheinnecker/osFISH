
pooling_strategy <- read_tsv("/g/schwab/rheinnec/projects/osFISH/dino_pooling_strategy - Sheet1.tsv")


primers <- tibble(
  
  primer_id=c(25, 27, 28),
  pool_id=c("p1", "p2", "p3"),
  primer_seq=c("CCAATAATA",
               "CATCATCAT",
               "CAACTTAAC")
  
)



required_probes <- pooling_strategy %>%
  select(species=`3pools_order`, p1, p2, p3) %>%
  pivot_longer(cols = c("p1", "p2", "p3"), names_to = "pool_id", values_to = "identifier") %>%
  filter(!is.na(identifier))




final_probes <- required_probes %>%
  left_join(anno_species) %>%
  select(all_of(names(required_probes)), revCom, Tm) %>%
  left_join(primers) %>%
  mutate(
    final_probe_seq=paste0(revCom, "TTT", primer_seq)
  )

write_tsv(final_probes, file="/g/schwab/Marco/projects/osFISH/probe_design/probes_selected_order001_05032025.tsv")
