library(tidyverse)
all_probes_files <- list.files(path="/g/schwab/Marco/projects/osFISH/species_seq", 
                               pattern="_probes_log.tsv", 
                               recursive = T,
                               full.names = T)


full_probes <- lapply(all_probes_files, function(file){
  
  read_tsv(file) %>%
    mutate(
      file=file
    )
  
}) %>%
  bind_rows()


anno_species <- full_probes %>%
  mutate(
    species=basename(dirname(file))
  )




anno_species %>%
  #group_by() %>%
  #filter(mismatch1_abs<1000) %>%
  group_by(species) %>%
  arrange(
    desc(mismatch1_abs),
    desc(class)
  ) %>%
  mutate(
    rank_mismatch1_abs=seq(1,length(species))
  ) %>%
  arrange(
    #desc(mismatch1_abs),
    desc(class)
  ) %>%
  mutate(
    rank_class=seq(1,length(species))
  ) %>%
  View()







