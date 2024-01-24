library(tidyverse)

df1 <- read_csv('deBruijn_graphs/keto/results/f_7mers_uniq.csv') %>%
  mutate(family = 'F')
df2 <- read_csv('deBruijn_graphs/keto/results/k_7mers_uniq.csv') %>%
  mutate(family = 'K')


foldchange <- bind_rows(df1, df2) %>%
  pivot_wider(names_from='family', values_from = 'count') %>%
  replace(is.na(.), as.numeric("0")) %>%
  mutate(F_impute = F+1, K_impute = K+1) %>%
  mutate(F_prop = (F/67)*100, K_prop = (K/167)*100) %>% 
  mutate(fc = F_impute/K_impute,
         log2fc = log2(fc)) %>%
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% 
  arrange(desc(log2fc)) %>%
  mutate(motif = paste0(node1, substr(node2, nchar(node2), nchar(node2)))) %>%
  select(node1, node2, motif, everything(), -ids)

write_csv(foldchange, 'deBruijn_graphs/keto/results/k-vs-f_motif_enrichment_060623.csv')

gp_df <- bind_rows(f_df, k_df) %>%
  group_by(node1, node2) %>%
  summarize(.groups="keep",
            paste(unique(family), collapse = ', '),
            weight = sum(count)) %>%
  arrange(desc(weight))

write_csv(fk_df, 'deBruijn_graphs/keto/results/fk_7mers_uniq.csv')

library(tidygraph)
library(ggraph)  

palette_pretty <- c("#009E24", "#0072B2","#E69F00", "#FF0000", "#979797", "#5530AA", "#1E1E1E")

# reproducibility seed
set.seed(13)

# create tidygraph object
tg <- as_tbl_graph(fk_df, directed = TRUE)

# plot graph using ggraph
p <- ggraph(tg, layout = "kk") + 
  geom_edge_link(aes(color = family, width = weight),
                 #aes(color = ProteinCount, alpha = ProteinCount),
                 #arrow = arrow(length = unit(0.05, "inches")),
                 #edge_width = 2,
                 linejoin = "round", 
                 lineend = "round",
  ) +
  scale_edge_width_continuous(range = c(1, 3)) +
  #scale_alpha_continuous(range = c(0.2,1)) +
  #xlim(-5, 5) +
  #ylim(-5, 5) +
  theme_void() +
  theme(legend.position="none") +
  #scale_edge_colour_gradient(low = "#132B43", high = "#56B1F7") +
  scale_edge_colour_viridis(direction = -1, begin = 0, end = 1, option = "rocket")
p