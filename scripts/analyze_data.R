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

oligo_1_files <- dir('oligos/results/oligo_1/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_2_files <- dir('oligos/results/oligo_2/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_3_files <- dir('oligos/results/oligo_3/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
strep_files <- dir('oligos/results/streptavidin/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)

strep_data <- combine_data(strep_files, fdr_cutoff = cutoff)
oligo_1_data <- combine_data(oligo_1_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))
oligo_2_data <- combine_data(oligo_2_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))
oligo_3_data <- combine_data(oligo_3_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label))


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

# plot log2fc for each label for each k
oligo_1_data %>%
  ggplot(aes(x = as.factor(kmer_size), y = `3w`, fill = label)) +
  geom_boxplot() +
  scale_fill_manual(values=pal1)

oligo_2_data %>%
  ggplot(aes(x = as.factor(kmer_size), y = `3w`, fill = label)) +
  geom_boxplot() +
  scale_fill_manual(values=pal1)

oligo_3_data %>%
  ggplot(aes(x = as.factor(kmer_size), y = `3w`, fill = label)) +
  geom_boxplot() +
  scale_fill_manual(values=pal1)

# what is unique to each data set?
#install.packages("ggVennDiagram")
library(ggVennDiagram)
k = 10
oligo_1_k_subset <- oligo_1_data %>% filter(kmer_size == k)
oligo_2_k_subset <- oligo_2_data %>% filter(kmer_size == k)
oligo_3_k_subset <- oligo_3_data %>% filter(kmer_size == k)
strep_k_subset <- strep_data %>% filter(kmer_size == k)
x <- list(oligo_1 = oligo_1_k_subset$edge, oligo_2 = oligo_2_k_subset$edge, 
          oligo_3 = oligo_3_k_subset$edge, streptavidin = strep_k_subset$edge)
ggVennDiagram(x)

ggVennDiagram(x[1:3], color = "black", lwd = 0.8, lty = 1, label = "count") +
  scale_fill_viridis_c(direction = -1, option = "A")

# compute PCA + plot
