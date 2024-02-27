library(tidyverse)
library(patchwork)
library(argparse)
theme_set(theme_bw(base_size = 12))
pal <- c("#E64B35", "#3578E6", "#4BD2AC", 
         "#665899", "#E6A435","#8A8A95", "#353540")

parser <- ArgumentParser()
parser$add_argument("-k", "--kmer_size", dest="kmer_size", required=TRUE,
                    help="Size of k for which to generate enrichment plots.")
args <- parser$parse_args()

# -----------------------------------
# functions
# -----------------------------------

combine_data <- function(data_dir, data_pattern){
  
  files <- dir(data_dir, pattern = paste0("*", data_pattern))
  
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

plot_log2fc <- function(df, k){
  
  p1 <- df %>% 
    filter(analyte != "streptavidin") %>% 
    ggplot(aes(y = log2fc, x = exp, fill = label)) +
    geom_boxplot(alpha = 0.75) +
    scale_fill_manual(values = pal) +
    theme(axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) +
    facet_wrap(~analyte, nrow = 1, scales = "free_x")
  p1
  
  out_prefix1 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_oligo_log2fc")
  h = 3; w = 7
  p1 %>% ggsave(paste0(out_prefix1, ".png"), ., device = "png",
                width = w, height = h, units = "in")
  p1 %>% ggsave(paste0(out_prefix1, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
  p2 <- df %>% filter(analyte == "streptavidin") %>%
    ggplot(aes(y = log2fc, x = exp, fill = label)) +
    geom_boxplot(alpha = 0.75) +
    scale_fill_manual(values = pal) +
    theme(axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) +
    facet_wrap(~analyte, nrow = 1, scales = "free_x")
  p2
  
  out_prefix2 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_strep_log2fc")
  h = 3; w = 3
  p2 %>% ggsave(paste0(out_prefix2, ".png"), ., device = "png",
                width = w, height = h, units = "in")
  p2 %>% ggsave(paste0(out_prefix2, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
}

plot_apex_zscore <- function(df, k){
  
  p1 <- df %>% 
    filter(analyte != "streptavidin") %>% 
    ggplot(aes(y = apex_zscore, x = exp, fill = label)) +
    geom_boxplot(alpha = 0.75) +
    scale_fill_manual(values = pal) +
    theme(axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) +
    facet_wrap(~analyte, nrow = 1, scales = "free_x")
  p1
  
  out_prefix1 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_oligo_zscore")
  h = 3; w = 7
  p1 %>% ggsave(paste0(out_prefix1, ".png"), ., device = "png",
                width = w, height = h, units = "in")
  p1 %>% ggsave(paste0(out_prefix1, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
  p2 <- df %>% filter(analyte == "streptavidin") %>%
    ggplot(aes(y = apex_zscore, x = exp, fill = label)) +
    geom_boxplot(alpha = 0.75) +
    scale_fill_manual(values = pal) +
    theme(axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) +
    facet_wrap(~analyte, nrow = 1, scales = "free_x")
  p2
  
  out_prefix2 = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/figures/enrich/", k, "mer_strep_zscore")
  h = 3; w = 3
  p2 %>% ggsave(paste0(out_prefix2, ".png"), ., device = "png",
                width = w, height = h, units = "in")
  p2 %>% ggsave(paste0(out_prefix2, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")
  
}


# -----------------------------------
# main
# -----------------------------------

k = args$kmer_size

kmer = paste0(k, "mer")
kmer_dir = paste0("/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/results/", kmer, "/enrichment_all/")
df <- combine_data(data_dir = kmer_dir,
                   data_pattern = "_enrich.csv") %>%
  filter(label != "library match") %>%
  mutate(analyte_exp = paste(analyte, exp, sep="_"))

plot_log2fc(df, k)
plot_apex_zscore(df, k)

