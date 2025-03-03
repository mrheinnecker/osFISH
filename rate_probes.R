library(tidyverse)

## if i did seperate fro each lengthm, use this blocm to combine files:

library(tidyverse)

PATT <- "_probes_log.tsv"
PATT <- "_probes_log_filtered_best4.tsv"
all_probes_files <- list.files(path="/scratch/rheinnec/osFISH/species/", 
                               pattern=PATT, 
                               recursive = T,
                               full.names = T)



lapply(unique(basename(dirname(all_probes_files))), function(SPEC){
  
  print(SPEC)
  
  comb_file <- lapply(all_probes_files[which(str_detect(all_probes_files, SPEC))], read_tsv) %>%
    bind_rows() %>%
    mutate(identifier=paste0(identifier, length))
  
  write_tsv(comb_file,
            file=file.path("/g/schwab/rheinnec/projects/osFISH/species_seq_40bp/", SPEC, paste0(SPEC, PATT)))
  
})

### from here probe analysis







all_probes_files <- list.files(path="/g/schwab/Marco/projects/osFISH/species_seq_40bp", 
                               pattern="_probes_log.tsv", 
                               recursive = T,
                               full.names = T)

all_best4_files <- list.files(path="/g/schwab/Marco/projects/osFISH/species_seq_40bp", 
                               pattern="best4.tsv", 
                               recursive = T,
                               full.names = T)


full_probes <- lapply(all_probes_files, function(file){
  
  read_tsv(file) %>%
    mutate(
      file=file
    )
  
}) %>%
  bind_rows()%>%
  mutate(
    species=basename(dirname(file))
  )

full_best4 <- lapply(all_best4_files, function(file){
  
  read_tsv(file, col_types = cols(.default = "c")) %>%
    mutate(
      file=file
    )
  
}) %>%
  bind_rows()%>%
  mutate(
    species=basename(dirname(file))
  )


anno_species <- full_probes %>%
  left_join(
    full_best4 %>% select(identifier, species) %>% mutate(best4=T)
  )


anno_species %>%
  filter(hitsR_abs<50) %>%
  group_by(species) %>%
  tally()

ranked <- anno_species %>%
#  filter(best4) %>%
  #group_by() %>%
  #filter(mismatch1_abs<1000) %>%
  group_by(species) %>%
  arrange(
    hitsR_abs,
    mismatch1_abs
    #desc(class)
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
  mutate(
    rank_sum=rank_mismatch1_abs+rank_class
  ) %>%
  arrange(rank_sum) %>%
  mutate(
    final_rank=seq(1,length(species))
  ) 
  #View()



check <- ranked %>%
  select(species, identifier, final_rank, best4, hitsT_abs, hitsR_abs, mismatch1_abs, mismatch2_abs, class, Tm,revCom) %>%
  filter(final_rank<5|best4) %>%
  arrange(species)


write_tsv(check, file="/g/schwab/Marco/projects/osFISH/ov_table_probes.tsv")





