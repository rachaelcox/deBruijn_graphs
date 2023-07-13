library(tidyverse)

g_df <- read_csv('deBruijn_graphs/keto/results/g_7mers_uniq.csv') %>%
  mutate(family = 'G')
p_df <- read_csv('deBruijn_graphs/keto/results/p_7mers_uniq.csv') %>%
  mutate(family = 'P')


foldchange <- bind_rows(g_df, p_df) %>%
  pivot_wider(names_from='family', values_from = 'count') %>%
  replace(is.na(.), as.numeric("0")) %>%
  mutate(G_impute = G+1, P_impute = P+1) %>%
  mutate(G_prop = (G/42)*100, P_prop = (P/40)*100) %>% 
  mutate(fc = G_impute/P_impute,
         log2fc = log2(fc)) %>%
  mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% 
  arrange(desc(log2fc))

write_csv(foldchange, 'deBruijn_graphs/keto/results/p-vs-g_motif_enrichment_060223.csv')

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