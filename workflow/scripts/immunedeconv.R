log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## ImmuneDeconv
library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

## Variables
cibersort_path <- snakemake@params[["cibersort_path"]] #'/cluster/projects/mcgahalab/ref/immunedeconv/cibersort'
mcp_path <- snakemake@params[["mcp_path"]] #'/cluster/projects/mcgahalab/ref/immunedeconv/mcp_counter'
#pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/INSPIRE/'
#setwd(pdir)

## Gene mapping
if(snakemake@params[["species"]] == 'homo_sapiens'){
  library("org.Hs.eg.db")
  genome <- org.Hs.eg.db
} else if(snakemake@params[["species"]] == 'mus_musculus'){
  library("org.Mm.eg.db")
  genome <- org.Mm.eg.db
}
txby <- keys(genome, 'ENSEMBL')
gene_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")


## Load in TPM matrix and format it with gene symbol rownames with no duplicates
tpm <- read.table(snakemake@input[[1]], header=T,
                  stringsAsFactors = F, check.names = F)
symbols <- gene_ids[tpm$gene]
na_idx  <- which(is.na(symbols))
dup_idx <- which(duplicated(symbols))
tpm_filt <- tpm[-unique(sort(c(na_idx, dup_idx))), ]
rownames(tpm_filt) <- symbols[-unique(sort(c(na_idx, dup_idx)))]
tpm_filt <- tpm_filt[,-1, drop=F]

ggplotit <- function(res){
  cols <- rainbow(n=nrow(res), s = 0.5, v=0.6)
  cols <- sample(cols, length(cols), replace = F)
  
  res %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_manual(values=cols) +
    scale_x_discrete(limits = rev(levels(res)))
}

## Deconvolution:
# Quantiseq
quantiseq <- deconvolute(tpm_filt, "quantiseq")
gg_qs <- ggplotit(quantiseq)

# MCP-counter
mcp_probesets <- read.table(file.path(mcp_path, "probesets.txt"),
                            sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
mcp_genes <- read.table(file.path(mcp_path, "genes.txt"),
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, 
                        colClasses = "character", check.names = FALSE)
mcp = deconvolute(tpm_filt, "mcp_counter",  probesets=mcp_probesets, genes=mcp_genes)
mcp_rescale <- mcp
mcp_rescale[,-1] <- apply(mcp[,-1], 2, function(i) i/sum(i))
gg_mcp <- ggplotit(mcp)

# CIBERSORT
set_cibersort_binary(file.path(cibersort_path, "CIBERSORT.R"))
set_cibersort_mat(file.path(cibersort_path, "LM22.txt"))
cibersort <- deconvolute(tpm_filt, "cibersort")
gg_cs <- ggplotit(cibersort)

## Merge the deconvs
qs_map = map_result_to_celltypes(quantiseq, quantiseq$cell_type, "cibersort") %>%
  t %>% melt
cs_map = map_result_to_celltypes(cibersort, cibersort$cell_type, "cibersort") %>% 
  t %>% melt
mcp_map = map_result_to_celltypes(mcp_rescale, mcp_rescale$cell_type, "cibersort") %>% 
  t %>% melt
maps <- list("quantiseq"=qs_map, "cibersort"=cs_map, "mcp"=mcp_map)
deconv_melt <- Reduce(function(x,y) merge(x, y, by=c('Var1', 'Var2'), all=T), maps)
colnames(deconv_melt)[-c(1:2)] <- names(maps)
deconv_melt <- melt(deconv_melt)
colnames(deconv_melt) <- c('sample', 'cell_type', 'method', 'score')

gg_mrg <- ggplot(deconv_melt, aes(x=sample, y=score, fill=method)) +
  facet_wrap(vars(cell_type), scales='free', ncol=3) +
  geom_bar(stat='identity', position='dodge') +
  ylim(0,1) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        strip.text.x = element_text(size = 8))

## Visualizations
pdf(snakemake@output[["aggregate"]], width=13, height = 13)
gg_mrg
dev.off()
pdf(snakemake@output[["mcp"]])
gg_mcp
dev.off()
pdf(snakemake@output[["quantiseq"]])
gg_qs
dev.off()
pdf(snakemake@output[["cibersort"]], width = 10)
gg_cs
dev.off()

## Data
saveRDS(maps, snakemake@output[["rds"]])
