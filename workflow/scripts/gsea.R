log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("org.Hs.eg.db")
library("KEGG.db")

library("GO.db")
library("Category")
library("GOstats")

library("clusterProfiler")
library("enrichplot")
library("ggplot2")


# Create a reference map of ENSEMBL to SYMBOL
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='ENTREZID',
                   keytype='ENSEMBL', multiVals="first")

## Read gene expr and DEG
res <- read.table(snakemake@input[[1]], sep="\t",
                  header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(res) <- res$gene

## Select the significant genes and filters
resSigind = res[ which(res$padj < as.numeric(as.character(snakemake@params[["maxp"]])) &
                         res$log2FoldChange > 0 &
                         res$baseMean > as.numeric(as.character(snakemake@params[["minbase"]]))), ]
resSigrep = res[ which(res$padj < as.numeric(as.character(snakemake@params[["maxp"]])) &
                         res$log2FoldChange < 0 &
                         res$baseMean > as.numeric(as.character(snakemake@params[["minbase"]]))), ]
resSig = rbind(resSigind, resSigrep)
resFilt <- res[res$baseMean > as.numeric(as.character(snakemake@params[["minbase"]])),]

## GO ontology enrichment analysis
print("GO Enrichment Analysis")
params=new("GOHyperGParams",
           geneIds=unique(na.omit(gene_ids[rownames(resSig)])),
           universeGeneIds=unique(na.omit(gene_ids[resFilt$gene])),
           annotation="org.Hs.eg.db",
           ontology="BP",
           pvalueCutoff=0.001,
           conditional=TRUE,
           testDirection="over")
overRepresented=hyperGTest(params)
go_summ <- summary(overRepresented)[,c(1,2,5,6,7)]

write.table(go_summ, file=snakemake@output[["go"]],
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

## GSEA
print("GSEA on GO terms")
gene_list <- resFilt$log2FoldChange
names(gene_list) <- resFilt$gene
gene_list<-sort(na.omit(gene_list), decreasing = TRUE)

gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "ENSEMBL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "fdr")

pdf(snakemake@output[["gseago_pdf"]], height = 20, width = 20)
print(emapplot(pairwise_termsim(gse), showCategory = 50))
print(ridgeplot(gse) + labs(x = "enrichment distribution"))
#gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
dev.off()

write.table(gse@result[,1:11], file=snakemake@output[["gseago_table"]],
            col.names = TRUE, row.names = TRUE, quote = FALSE)



## KEGG GSEA
print("GSEA on KEGG terms")
kegg_gene_list <- resFilt$log2FoldChange
names(kegg_gene_list) <- gene_ids[resFilt$gene]
kegg_gene_list <- kegg_gene_list[-which(is.na(names(kegg_gene_list)))]
kegg_gene_list<-sort(na.omit(kegg_gene_list), decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
             organism     = 'hsa',
             nPerm        = 10000,
             minGSSize    = 3,
             maxGSSize    = 800,
             pvalueCutoff = 0.05,
             pAdjustMethod = "fdr",
             keyType       = "ncbi-geneid",
             use_internal_data=TRUE)

pdf(snakemake@output[["gseakegg_pdf"]], height = 20, width = 20)
print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") +
facet_grid(.~.sign))
print(emapplot(pairwise_termsim(kk2), showCategory = 50))
print(ridgeplot(kk2) + labs(x = "enrichment distribution"))
dev.off()

write.table(kk2@result[,1:11], file=snakemake@output[["gseakegg_table"]],
          col.names = TRUE, row.names = TRUE, quote = FALSE)
