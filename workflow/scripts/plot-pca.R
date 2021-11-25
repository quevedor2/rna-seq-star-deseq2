log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggplot2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
tcnts <- as.data.frame(t(assay(counts)))
pca <- prcomp(tcnts, scale=F)
pca_x <- as.data.frame(cbind(pca$x[,c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')],
                       "condition"=as.character(counts[[snakemake@params[["pca_labels"]]]])))
for(id in paste0("PC", c(1:6))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
pca_x$Name <- rownames(pca_x)

pdf(snakemake@output[[1]])
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
ggplot(data=pca_x, aes(x=PC1, y=PC2, color=condition, label=Name)) +
      geom_point() + geom_text(aes(label=Name),hjust=0, vjust=0)
ggplot(data=pca_x, aes(x=PC3, y=PC4, color=condition, label=Name)) +
      geom_point() + geom_text(aes(label=Name),hjust=0, vjust=0)
ggplot(data=pca_x, aes(x=PC5, y=PC6, color=condition, label=Name)) +
      geom_point() + geom_text(aes(label=Name),hjust=0, vjust=0)
dev.off()
