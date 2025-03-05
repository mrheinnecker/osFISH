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

#write_tsv(anno_species, file = "/g/schwab/Marco/projects/osFISH/probe_design/oligoNdesign_per_species/combined_probes_all_species.tsv")



anno_species %>%
  filter(hitsR_abs<50) %>%
  group_by(species) %>%
  tally()

ranked <- anno_species %>%
#  filter(best4) %>%
  #group_by() %>%
  #filter(mismatch1_abs<1000) %>%
  group_by(species) %>%
  mutate(diff_tm=abs(66-Tm)) %>%
  filter(diff_tm<5) %>%
  arrange(
    hitsR_abs,
    mismatch1_abs,
    class
    #desc(class)
  ) %>%
  mutate(
    rank_mismatch1_abs=seq(1,length(species))
  ) 


## load manually selected prio 1 probes
it1_probes_raw <- read_tsv("/g/schwab/rheinnec/projects/osFISH/probes_it1.tsv", col_names = c("species", "identifier")) %>%
  mutate(it1=T) %>%
  left_join(ranked) #%>%

it1_probes <- it1_probes_raw %>%
  select(species, identifier, it1, it1_start_position=start_position)


pot_secondary_probe <- ranked %>%
  left_join(it1_probes %>% select(-it1_start_position), by=c("species", "identifier")) %>%
  filter(is.na(it1)) %>%
  filter(between(Tm, min(it1_probes_raw$Tm), max(it1_probes_raw$Tm))) %>%
  left_join(it1_probes %>% select(-identifier, -it1), by="species") %>%
  mutate(
    pos_diff=abs(start_position-it1_start_position)
  ) %>%
  filter(pos_diff>50) %>%
  group_by(species)
  



pot_secondary_probe %>%
  group_by(species) %>%
  tally()



best20 <- ranked %>% filter(rank_mismatch1_abs<3)

ranked %>% filter(rank_mismatch1_abs<30) %>%
  group_by(species, region) %>%
  tally() %>%
  group_by(species) %>%
  tally() %>%
  arrange(n)
  #View()
  

#%>%
  arrange(
    #desc(mismatch1_abs),
    class
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


write_tsv(check, file="/g/schwab/Marco/projects/osFISH/ov_table_probes_new.tsv")

ranked %>%
  filter(species %in% c("Heterocapsa_rotundata", "Karenia_mikimotoi"))



