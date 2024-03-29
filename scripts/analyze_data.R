library(tidyverse)
library(patchwork)
library(fs)
theme_set(theme_bw(base_size = 12))
pal1 <- c("#E64B35", "#3578E6", "#00A087", "#5A35E3", "#3C5488",
          "#353540", "#F39B7F", "#4DBBD5", "#91D1C2", "#E6A435")

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

calc_pca <- function(df){
  
  # do pca
  pca <- df %>%
    select(-edge, -confect, -label) %>%
    mutate_if(is.character, as.numeric) %>%
    scale() %>%
    prcomp()
  
  # add exp information back into pca data
  pca_data <- data.frame(pca$x, edge = df$edge, 
                         label = df$label)
}

# ----------------------------
# main
# ----------------------------

# reading in data
cutoff = 0.10

strep_files <- dir('oligos/results/streptavidin/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_1_files <- dir('oligos/results/oligo_1/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_2_files <- dir('oligos/results/oligo_2/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_3_files <- dir('oligos/results/oligo_3/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)
oligo_b_files <- dir('oligos/results/oligo_b/labeled_nolib', pattern='*all_washes_enriched_nolib_labeled.csv', full.names=TRUE)


strep_data <- combine_data(strep_files, fdr_cutoff = cutoff) %>%
  mutate(exp = "streptavidin")

oligo_1_data <- combine_data(oligo_1_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label)) %>%
  mutate(exp = "oligo 1")

oligo_2_data <- combine_data(oligo_2_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label)) %>%
  mutate(exp = "oligo 2")

oligo_3_data <- combine_data(oligo_3_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label)) %>%
  mutate(exp = "oligo 3")

oligo_b_data <- combine_data(oligo_b_files, fdr_cutoff = cutoff) %>%
  mutate(label = case_when(edge %in% strep_data$edge ~ 'streptavidin bound',
                           label == 'no match' ~ 'everything else',
                           TRUE ~ label)) %>%
  mutate(exp = "oligo b")


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

oligo_b_data %>%
  ggplot(aes(x = as.factor(kmer_size), fill = label)) +
  geom_bar(position = position_dodge2(width = 1, preserve = "single")) +
  scale_y_log10() +
  scale_fill_manual(values=pal1) +
  labs(title = "oligo b")

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

oligo_b_data %>%
  ggplot(., aes(x = as.factor(kmer_size), y = `3w`, fill = label)) +
  geom_boxplot() +
  scale_fill_manual(values=pal1)

# what is unique to each data set?
#install.packages("ggVennDiagram")
library(ggVennDiagram)
k = 15
oligo_1_k_subset <- oligo_1_data %>% filter(kmer_size == k)
oligo_2_k_subset <- oligo_2_data %>% filter(kmer_size == k)
oligo_3_k_subset <- oligo_3_data %>% filter(kmer_size == k)
strep_k_subset <- strep_data %>% filter(kmer_size == k)

x <- list(oligo_1 = oligo_1_k_subset$edge, 
          oligo_2 = oligo_2_k_subset$edge, 
          oligo_3 = oligo_3_k_subset$edge, 
          streptavidin = strep_k_subset$edge)

ggVennDiagram(x, color = "black", lwd = 0.8, lty = 1, 
              label = "count") +
  scale_fill_viridis_c(direction = -1, option = "B")

ggVennDiagram(x[1:3], color = "black", lwd = 0.8, lty = 1, 
              label = "count") +
  scale_fill_viridis_c(direction = -1, option = "A")

# compute PCA + plot
slice_data <- function(df, top_n, k, filter=NA){
  sliced <- df %>%
    filter(kmer_size == k, label != filter) %>% 
    slice_max(order_by = `3w`, n = as.numeric(top_n)) %>% 
    select(edge, avg_expr, fdr, `0w`, `3w`, kmer_size, exp, label) %>%
    arrange(desc(`3w`))
  return(sliced)
}

# ------------------------------
# take top kmers for each size of k
# ------------------------------

slice_all_k <- function(df, top_n, filter=""){
  
  sliced <- df %>%
    filter(label != filter,
           kmer_size %in% c(11,13,15)) %>%
    group_by(kmer_size) %>% 
    slice_max(order_by = `3w`, n = as.numeric(top_n)) %>% 
    mutate(p = strsplit(edge, "")) %>% 
    unnest_wider(p, names_sep = "_") %>%
    mutate_each(funs(chartr("ATCG","1234",.))) %>%
    select(edge, avg_expr, fdr, `0w`, `3w`, kmer_size, exp, label,
           matches("^p_\\d*")) %>%
    ungroup()
  
  make_numeric <- sliced %>% select(matches("^p_\\d*")) %>% colnames()
  sliced <- sliced %>% mutate_at(vars(make_numeric), as.numeric)
  sliced[is.na(sliced)] <- 0
  
  return(sliced)
}

top_n = 100
oligo_1_slice <- slice_all_k(oligo_1_data, top_n, filter="streptavidin bound")
oligo_2_slice <- slice_all_k(oligo_2_data, top_n, filter="streptavidin bound")
oligo_3_slice <- slice_all_k(oligo_3_data, top_n, filter="streptavidin bound")
oligo_b_slice <- slice_all_k(oligo_b_data, top_n, filter="streptavidin bound")
strep_slice <- slice_all_k(strep_data, top_n)

all_data <- rbind(oligo_1_slice, oligo_2_slice,
                      oligo_3_slice, strep_slice) 

all_data_pca <- all_data %>%
  select(-edge, -label, -kmer_size, -exp) %>% 
  select(matches("p_\\d*")) %>% 
  mutate_if(is.character, as.numeric) %>% 
  scale() %>%
  prcomp()

# add exp information back into pca data
pca_labeled <- data.frame(all_data_pca$x, 
                          edge = all_data$edge, 
                          label = all_data$label,
                          exp = all_data$exp)

# plot PC1 and PC2
pca_labeled %>%
  ggplot(., aes(x = PC1, y = PC2, color = exp)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = pal1) +
  theme(legend.position = "top",
        legend.title = element_blank())

percent <- 100*all_data_pca$sdev^2 / sum(all_data_pca$sdev^2)
perc_data <- data.frame(percent = percent, PC = 1:length(percent))
ggplot(perc_data, aes(x = PC, y = percent)) + 
  geom_col() +
  ylab("% variance explained")

# ------------------------------
# PCA
# ------------------------------
top_n = 50
k = 15
oligo_1_slice <- slice_data(oligo_1_data, top_n, k, 
                            filter="streptavidin bound")
oligo_2_slice <- slice_data(oligo_2_data, top_n, k,
                            filter="streptavidin bound")
oligo_3_slice <- slice_data(oligo_3_data, top_n, k,
                            filter="streptavidin bound")
strep_slice <- slice_data(strep_data, top_n, k)
  

all_data <- rbind(oligo_1_slice, oligo_2_slice, oligo_3_slice,
                  strep_slice) %>%
  mutate(p = strsplit(edge, "")) %>% 
  unnest_wider(p, names_sep = "_") %>%
  mutate_each(funs(chartr("ATCG","1234",.)))

pca <- all_data %>%
  select(-edge, -label, -kmer_size, -exp) %>% 
  select(matches("p_\\d*")) %>% 
  mutate_if(is.character, as.numeric) %>% 
  scale() %>%
  prcomp()

# add exp information back into pca data
pca_labeled <- data.frame(pca$x, 
                          edge = all_data$edge, 
                          label = all_data$label,
                          exp = all_data$exp)

# plot PC1 and PC2
pca_labeled %>%
  filter(exp != "oligo b") %>% 
  ggplot(., aes(x = PC1, y = PC2, color = exp)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = pal1) +
  theme(legend.position = "top",
        legend.title = element_blank())

# plot % variance explained
percent <- 100*pca$sdev^2 / sum(pca$sdev^2)
perc_data <- data.frame(percent = percent, PC = 1:length(percent))
ggplot(perc_data, aes(x = PC, y = percent)) + 
  geom_col() +
  ylab("% variance explained")
