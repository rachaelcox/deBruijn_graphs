library(tidyverse)
library(viridis)

setwd("/stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/species_codes/")

crossref_0 = read_csv("species-short_domains_taxons_formatted.txt")
crossref_0

oma_longID = read_csv("species-full_taxons_formatted.txt")
oma_longID

crossreferences <- oma_longID %>%
  left_join(crossref_0, by=c("TaxonID")) %>% 
  mutate(RE_species=paste0(OMA_full,"[0-9]"))


palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
ggplot(crossreferences) +
  geom_bar(aes(x = Domain, fill = Domain)) +
  #scale_fill_viridis() +
  scale_fill_manual(values = palette_OkabeIto) +
  theme_classic()

OMA_eukarya <- crossreferences %>%
  filter(Domain=="E")

OMA_bacteria <- crossreferences %>%
  filter(Domain=="B")

OMA_archaea <- crossreferences %>%
  filter(Domain=="A")

write_csv(OMA_eukarya, "OMA_eukarya.txt")
write_csv(OMA_bacteria, "OMA_bacteria.txt")
write_csv(OMA_archaea, "OMA_archaea.txt")
