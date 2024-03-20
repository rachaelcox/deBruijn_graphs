library(tidyverse)
library(patchwork)
library(fs)
theme_set(theme_bw(base_size = 12))
pal1 <- c("#E64B35", "#3578E6", "#00A087", "#3C5488","#353540",
          "#F39B7F", "#4DBBD5", "#91D1C2", "#E6A435")
pal2 <- c("#5A35E3", "#B435E3", "#3578E6", "#73787C", "#FDE725", 
                     "#4BD2AC", "#00A201")

# ----------------------------
# functions
# ----------------------------
combine_data <- function(files, fdr_cutoff){
  
  dfs <- list()
  for (k in seq.int(5, 15, 1)){
    
    idx <- grep(paste0('_', k, 'mers_'), files)
    print(paste0('Reading data: ', files[idx]))
    df <- read_csv(files[idx]) %>%
      filter(fdr <= fdr_cutoff) %>%
      mutate(kmer_size = k)
    dfs[[k]] <- df
  }
  print('Combining data ...')
  data <- reduce(dfs, rbind)
  print('Done!')
  return(data)
}

# ----------------------------
# main
# ----------------------------

# reading in data
cutoff = 0.10

oligo_1_files <- dir('oligos/results/oligo_1/labeled', pattern='*all_washes_enriched_labeled.csv', full.names=TRUE)
oligo_2_files <- dir('oligos/results/oligo_2/labeled', pattern='*all_washes_enriched_labeled.csv', full.names=TRUE)
oligo_3_files <- dir('oligos/results/oligo_3/labeled', pattern='*all_washes_enriched_labeled.csv', full.names=TRUE)
strep_files <- dir('oligos/results/streptavidin/labeled', pattern='*all_washes_enriched_labeled.csv', full.names=TRUE)

oligo_1_data <- combine_data(oligo_1_files, fdr_cutoff = cutoff) #%>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin match',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))
oligo_2_data <- combine_data(oligo_2_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin match',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))
oligo_3_data <- combine_data(oligo_3_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin match',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))
strep_data <- combine_data(strep_files, fdr_cutoff = cutoff)


# plot k-mer count for each label
oligo_1_data %>%
  ggplot(aes(x = as.factor(kmer_size), fill = label)) +
  geom_bar(position = position_dodge2(width = 1, preserve = "single")) +
  scale_y_log10() +
  scale_fill_manual(values=pal1) +
  labs(title = "oligo 1")

oligo_2_data %>%
  ggplot(aes(x = as.factor(kmer_size), fill = label)) +
  geom_bar(position = position_dodge2(width = 1, preserve = "single")) +
  scale_y_log10() +
  scale_fill_manual(values=pal1) +
  labs(title = "oligo 2")

oligo_3_data %>%
  ggplot(aes(x = as.factor(kmer_size), fill = label)) +
  geom_bar(position = position_dodge2(width = 1, preserve = "single")) +
  scale_y_log10() +
  scale_fill_manual(values=pal1) +
  labs(title = "oligo 3")

strep_data %>%
  ggplot(aes(x = as.factor(kmer_size), fill = label)) +
  geom_bar(position = position_dodge2(width = 1, preserve = "single")) +
  scale_y_log10() +
  scale_fill_manual(values=pal1) +
  labs(title = "streptavidin")

# plot average expression for each label
oligo_1_data %>%
  #filter(label != "library match") %>% 
  ggplot(aes(x = avg_expr, fill = label)) +
  geom_density(alpha = 0.7) +
  #scale_fill_manual(values=c(pal1[1],pal1[2],pal1[4],pal1[5])) +
  scale_fill_manual(values=pal1) +
  facet_wrap(~kmer_size, ncol=3, scales = "free_y")

# ----------------------------
# archive/testing
# ----------------------------
fdr_cutoff = 0.10
dfo <- read_csv('oligos/results/oligo_3/labeled/oligo_3_7mers_10c_all_washes_enriched_labeled.csv') %>%
  filter(fdr <= fdr_cutoff)
dfs <- read_csv('oligos/results/streptavidin/labeled/streptavidin_7mers_10c_all_washes_enriched_labeled.csv') %>%
  filter(fdr <= fdr_cutoff)

targets <- dfo %>%
  filter(label == "complement oligo match")

overlap <- dfo %>%
  filter(edge %in% dfs$edge)
