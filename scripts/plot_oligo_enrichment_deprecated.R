library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 12))
pal <- c("#009E24", "#0072B2", "#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")

# -----------------------------------
# functions
# -----------------------------------

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

process_exps <- function(pattern, exp_prefix, k, exp_type){
  
  file_pattern = paste0(exp_prefix, 
                        pattern, k, "mer_enrichment.csv")  # arg
  files <- dir(path = paste0('deBruijn_graphs/oligo/results/', k, 'mer'),
               pattern = file_pattern,
               full.names = TRUE)
  bg_file = dir(path = paste0('deBruijn_graphs/oligo/results/', k, 'mer'),
                pattern = paste0(exp_prefix, 
                                 "_bg_", k, "mer_enrichment.csv"),
                full.names = TRUE)
  
  print(paste0("files for ", k, "mers:"))
  print(files)
  
  df <- fmt_data(files) %>%
    mutate(label = ifelse(label == "off target", "no match", label))%>%
    select(-zscore)
  print(df)
  
  df_bg <- read_csv(bg_file) %>%
    rename(bg_count = count) %>%
    select(-exp, -zscore) %>%
    mutate(label = ifelse(label == "off target", "no match", label))
  print(df_bg)
  
  df_joined <- left_join(df, df_bg, by=c("node1", "node2", "label")) %>%
    filter(label %in% c("oligo match", "no match")) %>% 
    mutate(type = exp_type)
  
  print(df_joined)
  return(df_joined)
  
}

calc_fc <- function(df){
  
  df_norm <- df %>% 
    mutate(count = replace_na(count, 0),
           bg_count = replace_na(bg_count, 0),
           exp_pseudo = count + 1,
           ctrl_pseudo = bg_count + 1) %>%
    group_by(exp) %>%
    mutate(F0ctrl_pseudo = ctrl_pseudo / sum(ctrl_pseudo)) %>%
    mutate(F0exp_pseudo = exp_pseudo / sum(exp_pseudo)) %>%
    mutate(fc = case_when(F0ctrl_pseudo > F0exp_pseudo ~ -F0ctrl_pseudo/F0exp_pseudo,
                          F0ctrl_pseudo < F0exp_pseudo ~ F0exp_pseudo/F0ctrl_pseudo)) %>%
    mutate(log2fc = log2((F0exp_pseudo/F0ctrl_pseudo))) %>%
    # z-score calculations
    mutate(F0_ctrl = bg_count / (sum(bg_count))) %>%
    mutate(F0_exp = count / (sum(count))) %>%
    mutate(F1 = (bg_count + count) / (sum(bg_count) + sum(count))) %>%
    mutate(zscore =
             (F0_exp - F0_ctrl)/
             sqrt((
               ((F1 * (1-F1)) / sum(count))) +
                 ((F1 * (1-F1)) / sum(bg_count))
             ))
  
  return(df_norm)
  
}

summarize_kmers <- function(exp_prefix){
  
  dfs <- list()
  for (i in 5:15){
    
    print(i)
    
    data_fc <- process_exps("_\\dw_", exp_prefix, i, exp_type="direct") %>%
      calc_fc()
    
    data_comp_fc <- process_exps("_\\dw_comp_", exp_prefix, i, exp_type="comp") %>%
      calc_fc()
    
    df_comb <- bind_rows(data_fc, data_comp_fc)
    print(df_comb)
    
    df_all <- df_comb %>% 
      group_by(exp, label, type) %>%
      summarize(sum_nkc = sum(F0exp_pseudo),
                mean_nkc = mean(F0exp_pseudo),
                sum_nkc_bg = sum(F0ctrl_pseudo),
                mean_nkc_bg = mean(F0ctrl_pseudo),
                mean_nkc_log2fc = log2((mean_nkc/mean_nkc_bg))) %>%
      mutate(k = i)
    
    outname = paste0("deBruijn_graphs/oligo/results/summarized/",  
                     exp_prefix, "_", i, "mers_norm_summarized.csv")
    
    data.table::fwrite(df_comb, outname)
    
    dfs[[i]] <- df_all
    
    rm(data_fc)
    rm(data_comp_fc)
    gc()
    
  }
  
  df_final <- bind_rows(dfs)
  print(df_final)
  
  df_final <- df_final %>%
    mutate(exp_label = paste(exp, label, sep=", "),
           label_type = paste(label, type, sep=", "),
           exp_type = paste(exp, type, sep=", "),
    )
  
  return(df_final)
  
}

# -----------------------------------
# results for 5-15mers
# -----------------------------------

plot_log2fc_xkmer <- function(df_final, exp_prefix){
  
  blues <- rev(c("#6B0077","#73579B","#828BBC"))
  reds <- rev(c("#7D0112", "#A14C26", "#C17F45"))
  pal_br <- c(blues, reds)
  
  x <- unique(df_final$exp_type)
  exp_ord <- c(x[1], x[3], x[5], x[2], x[4], x[6])
  df_final$exp_type <- factor(df_final$exp_type, levels = exp_ord)
  
  x <- unique(df_final$exp_label)
  exp_ord <- c(x[1], x[3], x[5], x[2], x[4], x[6])
  df_final$exp_label <- factor(df_final$exp_label, levels = exp_ord)
  
  p <- df_final %>%
    filter(type == "comp") %>% 
    ggplot(aes(x = k, y = mean_nkc_log2fc, color = exp_label)) +
    geom_point() +
    geom_line(aes(group = exp_label)) +
    scale_color_manual(values = pal_br) +
    scale_x_discrete(limits=c(6, 8, 10, 12, 14)) +
    #facet_wrap(~type) +
    labs(color = "k-mer label",
         y = "log2fc (mean)")
  p
  
  plot_outname = paste0("deBruijn_graphs/oligo/figures/",  
                        exp_prefix, "_log2fc")
  p %>% ggsave(paste0(plot_outname, ".png"), ., device = "png",
                width = 5, height = 3, units = "in")
  p %>% ggsave(paste0(plot_outname, ".pdf"), ., device = "pdf",
                width = 5, height = 3, units = "in")
  
  return(p)
}

o1 <- summarize_kmers(exp_prefix = "oligo_1")
o2 <- summarize_kmers(exp_prefix = "oligo_2")
o3 <- summarize_kmers(exp_prefix = "oligo_3")
sa <- summarize_kmers(exp_prefix = "streptavidin")

p1 <- plot_log2fc_xkmer(o1, "oligo_1")
p2 <- plot_log2fc_xkmer(o2, "oligo_2")
p3 <- plot_log2fc_xkmer(o3, "oligo_3")
psa <- plot_log2fc_xkmer(sa, "streptavidin")
psa

# make panel from all oligo results
p1_lab <- p1 + labs(y = "mean log2fc")
p2_noy <- p2 + theme(axis.title.y = element_blank())
p3_noy <- p3 + theme(axis.title.y = element_blank())

log2fc_panel <- p1_lab + p2_noy + p3_noy + plot_layout(guides = 'collect')

panel_outname = paste0("deBruijn_graphs/oligo/figures/oligo_log2fc_panel")
log2fc_panel %>% ggsave(paste0(panel_outname, ".png"), ., device = "png",
              width = 9.7, height = 2.75, units = "in")
log2fc_panel %>% ggsave(paste0(panel_outname, ".pdf"), ., device = "pdf",
              width = 9.7, height = 2.75, units = "in")

# -------------------------------------------

# -------------------------------------------
# density plots
# -------------------------------------------

exp_prefix = "streptavidin"
log_cutoff = 3
df <- data.table::fread(paste0("deBruijn_graphs/oligo/results/10mer/enrichment_all/", exp_prefix, ""))

fc_pval <- fc %>%
  mutate(pval = pnorm(zscore, lower.tail = FALSE)) %>%
  mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval)))

pval_cutoff <- fc_pval %>%
  filter(pval <= 0.05) %>%
  summarize(min = min(log2fc))

fdr_cutoff <- fc_pval %>%
  filter(fdr_bh <= 1e-10) %>%
  summarize(min = min(log2fc))

dp1 <- ggplot(fc_pval, aes(x = log2fc, fill = label)) +
  geom_histogram(alpha = 0.75, bins = 250) + scale_y_log10() +
  #geom_density(alpha=0.75) +
  scale_fill_manual(values = c("#73579B", "#A14C26")) +
  #geom_vline(xintercept = pval_cutoff[,1], linetype = "dashed", alpha = 0.75) +
  #annotate("text", x=0.72*pval_cutoff[,1], y=0.75, label="p <= 0.05", angle=90, size=3) +
  #geom_vline(xintercept = fdr_cutoff[,1], linetype = "dashed", alpha = 0.75) +
  #annotate("text", x=10*fdr_cutoff[,1], y=0.7, label="FDR=1x10^-10", angle=90, size=2) +
  #geom_vline(xintercept = log_cutoff, linetype = "dotted", alpha = 0.75) +
  #annotate("text", x=1.1*log_cutoff, y=0.7, label="log2fc=3", angle=90, size=2) +
  theme(legend.position = "top",
        #legend.key.size = ,
        legend.title = element_blank()) +
  facet_wrap(~exp, ncol = 1) +
  labs(y = "# of 10mers")
dp1

plot_outname = paste0("deBruijn_graphs/oligo/figures/",  
                      exp_prefix, "_10mer_log2fc_density")

dp1 %>% ggsave(paste0(plot_outname, ".png"), ., device = "png",
              width = 5, height = 5, units = "in")
dp1 %>% ggsave(paste0(plot_outname, ".pdf"), ., device = "pdf",
              width = 5, height = 5, units = "in")


# -------------------------------------------
# volcano plots
# -------------------------------------------

exp_prefix = "oligo_2"
fc <- data.table::fread(paste0("deBruijn_graphs/oligo/results/summarized/", exp_prefix, "_10mers_norm_summarized.csv")) %>%
  filter(type == "comp") %>%
         #exp == "3w"
  mutate(pval = pnorm(zscore, lower.tail = FALSE)) %>%
  mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval)))

# thres_outname = paste0("deBruijn_graphs/oligo/results/thresholded/",
#                       exp_prefix, "_10mer_graph_fdr1e-10_3log2fc.csv")
# fc_pval <- fc %>%
#   mutate(pval = pnorm(zscore, lower.tail = FALSE)) %>%
#   mutate(fdr_bh = p.adjust(pval, method = "BH", n = length(pval))) %>%
#   filter(fdr_bh <= 1e-10 & log2fc >= 3) %>%
#   select(node1, node2, log2fc) %>% 
#   write_csv(thres_outname)

plot_outname = paste0("deBruijn_graphs/oligo/figures/",  
                      exp_prefix, "_10mer_volcano")
int = -log10(1e-10)
vp <- ggplot(fc, aes(x = log2fc, y = -log10(fdr_bh), color = label)) +
  geom_point(size=1.25, alpha=0.8) +
  scale_color_manual(values = c("#73579B", "#A14C26")) +
  facet_wrap(~exp) +
  geom_vline(xintercept = 3, linetype = "dashed", 
             alpha = 0.75, size = 1.25) +
  geom_hline(yintercept = int, linetype = "dashed", 
             alpha = 0.75, linewidth = 1.25) +
  guides(color = guide_legend(reverse=TRUE))

vp %>% ggsave(paste0(plot_outname, ".png"), ., device = "png",
              width = 10, height = 3.5, units = "in")
