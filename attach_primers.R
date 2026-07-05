library(tidyverse)
pooling_strategy <- read_tsv("/g/schwab/marco/projects/osFISH/dino_pooling_strategy - Sheet1.tsv")


primers <- tibble(
  
  primer_id=c(25, 27, 28),
  pool_id=c("p1", "p2", "p3"),
  primer_seq=c("CCAATAATA",
               "CATCATCAT",
               "CAACTTAAC"),
  id_primer=c("p25", "p27", "p28")
  
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
    final_probe_seq=paste0(revCom, "TTT", primer_seq),
    name=paste(species, identifier, sep="_")
  )

write_tsv(final_probes, file="/g/schwab/Marco/projects/osFISH/probe_design/probes_selected_order001_05032025.tsv")



write_tsv(final_probes %>% select(name, final_probe_seq), 
          col_names = F,
          file="/g/schwab/Marco/projects/osFISH/probe_design/probes_selected_order001_05032025_idt_upload.tsv")





## for secondary order.... just attaching but as doublets
## load all probes
comb_probes <- read_tsv("/g/schwab/marco/projects/osFISH/probe_design/combined_probes_05032025.tsv")


rel_probes_01 <- c(
  ## pair AS - GI
  "Akashiwo_sanguinea_oligoN2_p28",
  "Gymnodinium_impudicum_oligoN2456_p27",
  ## pair PL - TH
  "Pyrocystis_lunula_oligoN717_p28",
  "Takayama_helix_oligoN3387_p27",
  ## pair AC - PM
  "Prorocentrum_micans_oligoN2271_p28",
  "Amphidinium_carterae_oligoN45337_p27",
  ## HR - Kmiki
  "Heterocapsa_rotundata_oligoN700_p27",
  "Karenia_mikimotoi_oligoN1323_p28",
  ## P - M
  "Protodinium_oligoN2836_p27",
  "Margalefidinium_oligoN61_p28"
)

rel_probes_02 <- c(
  ## pair AS - GI
  "Akashiwo_sanguinea_oligoN67_p25",
  "Akashiwo_sanguinea_oligoN17_p28",
  "Gymnodinium_impudicum_oligoN10_p25",
  # ## pair PL - TH
  "Pyrocystis_lunula_oligoN95_p28"
  # "Takayama_helix_oligoN3387_p27",
  # ## pair AC - PM
  # "Prorocentrum_micans_oligoN2271_p28",
  # "Amphidinium_carterae_oligoN45337_p27",
  # ## HR - Kmiki
  # "Heterocapsa_rotundata_oligoN700_p27",
  # "Karenia_mikimotoi_oligoN1323_p28",
  # ## P - M
  # "Protodinium_oligoN2836_p27",
  # "Margalefidinium_oligoN61_p28"
)



to_order <- tibble(
  id_full=rel_probes_02
) %>%
  mutate(
    id=str_remove(rel_probes, "_p27$|_p28$|_p25$"),
    id_primer=str_extract(id_full, "p27$|p28$|p25$")
  )  %>%
  left_join(comb_probes%>%
  mutate(id=paste(species, identifier, sep="_")), by="id")




final_probe_sequences <- to_order %>%
  ## adding scrambled sequence... no hits in BLAST and PR2 database
  # bind_rows(tibble(id="scrambled_N1", revCom="CTACGGTGGGACGAGAAATCTTACACAATCTGTGCGAAGT", 
  #                  id_primer=c("p27", "p28"))) %>%
  select(id, revCom, id_primer)%>%
  left_join(primers, by="id_primer") %>%
  rowwise() %>%
  mutate(
    seq_full=
      paste0(
        ## target sequence
        revCom,
        ## linker sequence
        "TTTTT",
        "A",
        primer_seq,
        "A",
        primer_seq
      ),
    short_id=case_when(
      str_detect(id, "mikimotoi") ~ "Kmiki",
      TRUE ~ str_split(id, "_") %>% 
      lapply(str_sub,1,1) %>% unlist() %>% .[1:(length(.)-1)] %>% toupper() %>% paste(collapse="")
    ),
    order_id=paste(short_id, str_remove(str_extract(id, "N\\d+"), "N"), sep="") %>%
      paste0(id_primer)
  )




write_tsv(final_probe_sequences %>% select(order_id, seq_full),
          file="/g/schwab/marco/projects/osFISH/probe_design/probes_selected_order02_20251111.tsv")







## select probes of relevance
AS <- comb_probes %>%
  filter(species=="Akashiwo_sanguinea")

PL <- comb_probes %>%
  filter(species=="Pyrocystis_lunula")

AC <- comb_probes %>%
  filter(species=="Amphidinium_carterae")

M <- comb_probes %>%
  filter(species=="Margalefidinium_oligoN61")


GI <- comb_probes %>%
  filter(species=="Gymnodinium_impudicum")






