library(tidyverse)
library(patchwork)
library(argparse)
theme_set(theme_bw(base_size = 12))
pal <- c("#E64B35", "#3578E6", "#4BD2AC", "#E6A435",
         "#665899","#8A8A95", "#353540")

parser <- ArgumentParser()
parser$add_argument("-k", "--kmer_size", dest="kmer_size", required=TRUE,
                    help="Size of k for which to generate enrichment plots.")
parser$add_argument("-n", "--top_n", dest="top_n", required=TRUE,
                    help="Number of top k-mers to consider.")
args <- parser$parse_args()

# -----------------------------------
# functions
# -----------------------------------
combine_data <- function(data_dir, data_pattern){
  
  files <- dir(data_dir, pattern = paste0("*", data_pattern))
  print(files)
  
  df <- tibble(filename = files) %>% 
    mutate(file_contents = map(filename, ~ data.table::fread(
      file.path(data_dir, .)))
    ) %>% 
    unnest() %>%
    mutate(filename = str_remove(filename, data_pattern)) %>% 
    mutate(analyte = case_when(
      str_detect(filename, "oligo_b") ~ "oligo_b",
      str_detect(filename, "oligo_1") ~ "oligo_1",
      str_detect(filename, "oligo_2") ~ "oligo_2",
      str_detect(filename, "oligo_3") ~ "oligo_3",
      str_detect(filename, 
                 "streptavidin") ~ "streptavidin")
    )
  
}


plot_pca <- function(df, k){
  
  # do pca
  pca <- df %>%
    select(-filename,-label, -analyte, -node1, -node2,
           -matches("*pval*")) %>%
    mutate(exp = gsub("w", "", exp)) %>% 
    mutate_if(is.character, as.numeric) %>%
    #select(matches("p_*")) %>% 
    scale() %>%
    prcomp()
  
  # add exp information back into pca data
  pca_data <- data.frame(pca$x, analyte = df$analyte, exp = df$exp, 
                         label = df$label)
  
  # plot pca
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = analyte)) + 
    geom_point(size = 3, alpha = 0.55) +
    scale_color_manual(values = pal) +
    scale_x_log10() +
    theme(legend.position = "top") +
    guides(color = guide_legend(nrow = 2))
  pca_plot
  
  h = 4.5; w = 4.5
  out_prefix1 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_pca_plot")
  pca_plot %>% ggsave(paste0(out_prefix1, ".png"), ., device = "png",
                width = w, height = h, units = "in")
  pca_plot %>% ggsave(paste0(out_prefix1, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
  # % variance explain
  # calculate variance
  percent <- 100*pca$sdev^2 / sum(pca$sdev^2)
  perc_data <- data.frame(percent = percent, PC = 1:length(percent))
  perc_plot <- ggplot(perc_data, aes(x = PC, y = percent)) + 
    geom_col() +
    ylab("% variance explained")
  perc_plot
  
  h = 3; w = 4.5
  out_prefix2 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_pca_variance")
  perc_plot %>% ggsave(paste0(out_prefix2, ".png"), ., 
                       device = "png",
                      width = w, height = h, units = "in")
  perc_plot %>% ggsave(paste0(out_prefix2, ".pdf"), ., 
                       device = "pdf",
                      width = w, height = h, units = "in")
  
}


# -----------------------------------
# main
# -----------------------------------

k=args$kmer_size
top_n=args$top_n
#k=14
#top_n=100
kmer = paste0(k, "mer")
kmer_dir = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/results/", kmer, "/enrichment_all/")
df <- combine_data(data_dir = kmer_dir,
                   data_pattern = "_enrich_z_filter.csv") %>%
  filter(label != "library match", exp != "0w") %>%
  mutate(p = strsplit(node_edge, "")) %>% 
  unnest_wider(p, names_sep = "_") %>%
  mutate_each(funs(chartr("ATCG","1234",.))) %>%
  group_by(analyte, exp) %>% 
  slice_max(order_by = apex_zscore, n = as.numeric(top_n)) %>% 
  ungroup()

plot_pca(df, k)