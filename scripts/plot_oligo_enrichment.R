library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 12))
pal <- c("#009E24", "#0072B2","#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")

fmt_data <- function(files){
  
  dfs <- list()
  
  for (i in seq_along(files)){
    print(paste(i,files[[i]]))
    df <- read_csv(files[[i]])
    dfs[[i]] <- df
  }
  
  final_df <- bind_rows(dfs)
  
  return(final_df) 
  
}

exp_title = "Oligo-1"
exp_prefix = "oligo_1"
file_pattern = paste0(exp_prefix, 
                      "_\\dw_comp_10mer_enrichment.csv")
files <- dir(path = 'deBruijn_graphs/oligo/results/10mer',
             pattern = file_pattern,
             full.names = TRUE)

df <- fmt_data(files)
df <- df %>%
  mutate(label = ifelse(label == "off target", "no match", label))

p1 <- ggplot(df, aes(x = label, y = count, fill = exp)) +
  geom_boxplot() +
  scale_fill_manual(values = pal) +
  scale_y_log10() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        ) +
  labs(y = "# unique NGS k-mers (k=10)")
p1

# p1 %>% ggsave("deBruijn_graphs/oligo/figures/oligo_3_10mer_counts.png", ., device = "png", 
#        width = 6, height = 5, units = "in")


dff <- df %>%
  filter(label %in% c("oligo match", "no match"))

p2 <- ggplot(dff, aes(x = count, fill = label)) +
  geom_density_ridges(alpha=0.75, (aes(y = label))) +
  scale_fill_manual(values = pal[5:length(pal)]) +
  facet_wrap(~exp, ncol=1) +
  scale_x_log10() +
  geom_vline(xintercept = 25, alpha = 0.5, linetype = "dashed") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "# unique NGS k-mers (k=10)")
p2

panel <- p1 + p2 + plot_annotation(exp_title)
panel

panel_outname = paste0("deBruijn_graphs/oligo/figures/",
                 exp_prefix,
                 "_comp_10mer_panel.png")

panel %>% ggsave(panel_outname, ., device = "png",
                     width = 10.5, height = 5.5, units = "in")

out_prefix = paste0("deBruijn_graphs/oligo/results/",
                    exp_prefix)

filtered_outfile = paste0(out_prefix,
                          "_comp_10mers_noflanks_25r.csv")
washed_outfile = paste0(out_prefix,
                          "_comp_10mers_noflanks_25r_washed.csv")
                          
dff %>%
  select(node1, node2, count) %>%
  filter(count >= 25) %>% 
  write_csv(filtered_outfile)

dff %>%
  filter(exp %in% c("1w", "2w"),
         count >= 25) %>% 
  select(node1, node2, count) %>%
  write_csv(washed_outfile)