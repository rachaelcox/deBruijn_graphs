library(tidyverse)
library(patchwork)
library(scales)
theme_set(theme_bw(base_size = 12))
options(scipen=100000)
pal <- c("#009E24", "#0072B2", "#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")
pub_pal <- c("#E64B35", "#3578E6", "#4BD2AC", "#665899", "#E6A435","#8A8A95", "#353540")

# -----------------------------------
# functions
# -----------------------------------

combine_data <- function(data_dir, data_pattern){
  
  files <- dir(data_dir, pattern = paste0("*", data_pattern))
  
  df <- tibble(filename = files) %>% 
    mutate(file_contents = map(filename, ~ read_csv(
      file.path(data_dir, .)))
    ) %>%
    unnest() %>%
    mutate(filename = str_remove(filename, data_pattern)) %>% 
    mutate(exp = case_when(
      str_detect(filename, "oligo_b") ~ "oligo_b",
      str_detect(filename, "oligo_1") ~ "oligo_1",
      str_detect(filename, "oligo_2") ~ "oligo_2",
      str_detect(filename, "oligo_3") ~ "oligo_3",
      str_detect(filename, 
                 "streptavidin") ~ "streptavidin")
    )
  
}

plot_counts <- function(df, k){
  
  kmer = paste0(k, "mer")
  xlabel = paste("unique", kmer, "counts")
  
  p1 <- ggplot(df, aes(x = count, y = fct_rev(filename), fill = exp)) +
    geom_boxplot() +
    scale_fill_manual(values = pub_pal) +
    theme(axis.title.y = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) +
    labs(x = xlabel) +
    scale_x_log10(label=comma)
  
  out_prefix = paste0("deBruijn_graphs/oligo/figures/", kmer, "_counts")
  h = 6.5; w = 7.5
  p1 %>% ggsave(paste0(out_prefix, "_full.png"), ., device = "png",
                width = w, height = h, units = "in")
  p1 %>% ggsave(paste0(out_prefix, "_full.pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
}


# -----------------------------------
# main
# -----------------------------------

# for (k in 5:15){
#   
#   kmer = paste0(k, "mer")
#   kmer_dir = paste0("deBruijn_graphs/oligo/results/", kmer, "/")
#   kmer_pattern = paste0(kmer, "s_unique.csv")
#   df <- combine_data(data_dir = kmer_dir,
#                      data_pattern = kmer_pattern)
#   plot_counts(df, k=k)
# 
# }

k=15
kmer = paste0(k, "mer")
kmer_dir = paste0("deBruijn_graphs/oligo/results/", kmer, "/")
kmer_pattern = paste0(kmer, "s_unique.csv")
df <- combine_data(data_dir = kmer_dir,
                   data_pattern = kmer_pattern)
plot_counts(df, k=k)