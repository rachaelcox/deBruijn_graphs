library(limma)
library(edgeR)
library(jsonlite)
library(argparse)
library(tidyverse)
parser <- ArgumentParser()
parser$add_argument("-f", "--counts_file", dest="counts_file", required=TRUE,
                    help="File containing k-mer counts in wide format.")
parser$add_argument("-e", "--exp", dest="exp", required=TRUE,
                    help="Experiment or list of experiments to test (wash prefix; e.g., for 1 experiment: '0w','3w'; for multiple experiments: '0w,1w,3w'.")
parser$add_argument("-o", "--output_file", dest="output_file", required=TRUE,
                    help="Name of file to write results to.")
args <- parser$parse_args()

# ----------------------------------
# code for testing
# ----------------------------------
# counts_file <- "oligos/results/oligo_1/counts_combined/oligo_1_10mers_combined_10c_wide.csv"
# output_file <- "oligos/results/oligo_1/enrichment/oligo_1_10mers_10c_3w_enriched.csv"
# exp <- "3w"

# ----------------------------------
# main
# ----------------------------------
counts_file <- args$counts_file
output_file <- args$output_file
exp <- args$exp

# fileEncoding='UTF-8-BOM' should strip the BOM marker FEFF that some windows tools add
print("Reading in data ...")
x <- read.delim(counts_file,
              sep=",",
              check.names=FALSE,
              colClasses='character',
              na.strings=c(),
              skip=0,
              fileEncoding='UTF-8-BOM')

exp_list <- str_split(exp, ',')[[1]]
cols <- colnames(x)
design_cols <- c("bg", exp_list)
exp_cols <- cols[grepl(paste(exp_list, collapse = "|"), colnames(x))]
bg_cols <- cols[grep(paste0("^bg*"), cols)]
count_cols <- c(bg_cols, exp_cols)

print("Constructing design and contrast matrices ...")
design_list <- list()
contrast_list <- list()
for (i in design_cols){
  idx <- grep(i, count_cols)
  # design vector
  d <- c(rep(0, length(count_cols)))
  d[idx] <- 1
  design_list[[i]] <- d
  # contrast vector
  if (i != "bg"){
    idx <- grep(i, design_cols)
    c <- c(-1, rep(0, length(design_cols)-1))
    c[idx] <- 1
    contrast_list[[i]] <- c
  }
}

design <- do.call(cbind, design_list)
rownames(design) <- count_cols
print("Design matrix:")
print(design)

cont.matrix <- do.call(cbind, contrast_list)
rownames(cont.matrix) <- c("bg", exp_list)
print("Contrast matrix:")
print(cont.matrix)

export_cols <- c("edge", count_cols)

# now re-read the first header line
# workaround R problem that always has strip.white=T for read.table
colnames(x) <- scan(counts_file,
                    what="",
                    sep=",",
                    nlines=1,
                    strip.white=F,
                    quote = "\"",
                    skip="0",
                    fileEncoding='UTF-8-BOM'
)

# force numeric count columns
x[,count_cols] <- apply(x[,count_cols], 2, function(v) as.numeric(v)) 
counts <- x[, count_cols]

# keep rows based on string based filters of columns; rows must match all filters
x <- x %>% select("edge", all_of(exp_cols), all_of(bg_cols))

# ** count threshold args here
# ** TO DO: add these args
# for now, default behavior is no CPM-based filtering

filter_rows <- fromJSON('[]')
if (length(filter_rows)>0) {
  keepRows <- apply(apply(filter_rows, 1, function(r) grepl(r['regexp'], x[,r['column']], perl=T, ignore.case=T)), 1, all)
} else {
  keepRows <- rep(TRUE, nrow(x))
}

# keep only rows with this x_min minimum
x_min = 0
keepMin <- apply(counts, 1, max) >= x_min
# keep only rows with cpm above y_min in at least z_min samples
y_min = 0
z_min = 0
keepCpm <- rowSums(cpm(counts) > y_min) >= z_min

keep <- keepMin & keepCpm & keepRows
x <- x[keep,]
counts <- counts[keep,]

# TMM normalization
print("Applying TMM normalization ...")
nf <- calcNormFactors(counts)
y <- voom(counts, design, plot=FALSE,lib.size=colSums(counts)*nf)

print("Fitting model ...")
fit <- lmFit(y,design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# format results
out <- topTable(fit2, n=Inf, sort.by='none')

out2 <- cbind(fit2$coef,
              out[, c('P.Value','adj.P.Val','AveExpr')],
              x[, export_cols]) %>%
  rename(pval = P.Value, fdr = adj.P.Val, avg_expr = AveExpr) %>%
  select(edge, everything())

# ** FDR arg here
# ** TO DO: add this arg
# for now, default FDR = 10%
print("Calculating top confidence effectors ...")
fdr = 0.1
confect <- topconfects::limma_confects(fit2, coef=1, fdr=fdr)

# ** outfile args here
print("Writing results ...")
out2 <- cbind(out2, confect=confect$table$confect[order(confect$table$index)])

out2 %>% write_csv(output_file)

print("Done!")