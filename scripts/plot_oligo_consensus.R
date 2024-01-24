library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 12))
pal <- c("#009E24", "#0072B2","#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")

k = 10
conditions = c('fdr1e-10', 'fdr1e-10_3log2fc', 'fdr1e-10_sd', 'fdr1e-10_3log2fc_sd')
c = conditions[4]
file_pattern = paste0("*", k, "mer_top100_consensus_alignment_", c, ".csv")  # arg
files <- dir(path = 'deBruijn_graphs/oligo/results/reconstructions/',
             pattern = file_pattern,
             full.names = TRUE)

data <- files %>%
  map(read_csv) %>%
  reduce(rbind)
data

random_lev = mean(data$random_levenshtein_ratio)

ggplot(data, aes(x = rank, y = levenshtein_ratio, color = exp, group = exp)) +
  geom_line(size=1.5) +
  #geom_point(size=2) +
  scale_color_manual(values = pal) +
  labs(y='Levenshtein ratio\n(true oligo sequence vs reconstructed consensus)',
       x = 'Seed 10mer rank (log2fc)') +
  geom_hline(yintercept = random_lev, linetype="dotted", linewidth=1) +
  annotate("text", x=0, y=random_lev*0.97, label="Random", size=2) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "in")) -> p1
p1

p1 %>% ggsave(paste0("deBruijn_graphs/oligo/figures/",
                     "top100_lev_ratio_", c, ".png"),
              ., device = "png", width = 4.5, height = 4.5, units = "in")
# p1 %>% ggsave(paste0("deBruijn_graphs/oligo/figures/",
#                      "top100_lev_ratio_", c, ".pdf"),
#               ., device = "pdf", width = 4.5, height = 4.5, units = "in")


