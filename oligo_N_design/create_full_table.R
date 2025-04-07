### flora probes




flora_full_probe_table_best4 <- list.files("/scratch/rheinnec/flora_oligo/", pattern = "best4", full.names = T, recursive=T) %>%
  lapply(function(FILE){
    
    read_tsv(FILE, col_types = cols(.default = "c")) %>%
      mutate(target=basename(FILE) %>% str_replace("_probes.tsv", "")) %>%
      return()
    
  }) %>%
  bind_rows() %>%
  mutate(species=str_replace_all(target, "_probes_log_filtered_best4.tsv", "")) %>%
  select(-target)





flora_full_probe_table <- list.files("/scratch/rheinnec/flora_oligo", pattern = "log_filtered.tsv", full.names = T, recursive=T) %>%
  lapply(function(FILE){
    
    read_tsv(FILE, col_types = cols(.default = "c")) %>%
      mutate(target=basename(FILE) %>% str_replace("_probes.tsv", "")) %>%
      return()
    
  }) %>%
  bind_rows() %>%
  mutate(species=str_replace_all(target, "_probes_log_filtered.tsv", "")) %>%
  select(-target)




write_tsv(flora_full_probe_table, file="/g/schwab/Marco/flora_full_probe_tabl.tsv")







