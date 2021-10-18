log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

print(str(snakemake@params[['contrast']]))
contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
reslfc <- lfcShrink(dds, contrast=contrast, type="normal")
# sort by p-value
res <- res[order(res$padj),]
reslfc <- reslfc[order(reslfc$padj),]
# TODO explore IHW usage


# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
plotMA(reslfc, ylim=c(-2,2))
dev.off()

# Map ENSEMBL IDs to HUGO Symbols
# Create a reference map of ENSEMBL to SYMBOL
if(snakemake@params[["species"]] == 'homo_sapiens'){
  library("org.Hs.eg.db")
  genome <- org.Hs.eg.db
} else if(snakemake@params[["species"]] == 'mus_musculus'){
  library("org.Mm.eg.db")
  genome <- org.Mm.eg.db
} else {
  library("org.Hs.eg.db")
  genome <- org.Hs.eg.db
}
txby <- keys(genome, 'ENSEMBL')
gene_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

write.table(data.frame("gene"=rownames(res), "symbol"=gene_ids[rownames(res)], res),
            file=snakemake@output[["table"]], row.names=FALSE, sep='\t')
