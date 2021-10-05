log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("KEGG.db")
library("GO.db")
library("Category")
library("GOstats")

library("clusterProfiler")
library("enrichplot")
library("ggplot2")


# Create a reference map of ENSEMBL to SYMBOL
if(snakemake@params[['genome']] == 'mus_musculus'){
    library("org.Mm.eg.db")
    annotation <- "org.Mm.eg.db"
    genome <- org.Mm.eg.db
    org_code <- 'mmu'
} else if(snakemake@params[['genome']] == 'homo_sapiens'){
    library("org.Hs.eg.db")
    annotation <- "org.Hs.eg.db"
    genome <- org.Hs.eg.db
    org_code <- 'hsa'
} else {
    stop("Your config 'genome' should be set to either 'homo_sapiens' or 'mus_musculus'")
}
txby <- keys(genome, 'ENSEMBL')
gene_ids <- mapIds(genome, keys=txby, column='ENTREZID',
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
if(nrow(resSig) > 0 & nrow(resFilt) > 0){
  print("GO Enrichment Analysis")
  params=new("GOHyperGParams",
             geneIds=unique(na.omit(gene_ids[rownames(resSig)])),
             universeGeneIds=unique(na.omit(gene_ids[resFilt$gene])),
             annotation=annotation,
             ontology="BP",
             pvalueCutoff=0.001,
             conditional=TRUE,
             testDirection="over")
  overRepresented=hyperGTest(params)
  go_summ <- summary(overRepresented)[,c(1,2,5,6,7)]
} else {
  print("No significant genes to run: GO Enrichment Analysis")
  go_summ <- data.frame()
}
write.table(go_summ, file=snakemake@output[["go"]],
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

## GSEA
if(nrow(resFilt) > 0){
    print("GSEA on GO terms")
    gene_list <- resFilt$log2FoldChange
    names(gene_list) <- resFilt$gene
    gene_list<-sort(na.omit(gene_list), decreasing = TRUE)

    gse <- tryCatch({
      gseGO(geneList=gene_list,
                 ont ="ALL",
                 keyType = "ENSEMBL",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 OrgDb = genome,
                 pAdjustMethod = "fdr")
     }, error=function(e){NULL})
} else {
    print("No significant genes to run: GSEA on GO terms")
    gse <- NULL
}

pdf(snakemake@output[["gseago_pdf"]], height = 20, width = 20)
if(is.null(gse)){
  write.table(data.frame(), file=snakemake@output[["gseago_table"]],
              col.names = TRUE, row.names = TRUE, quote = FALSE)
} else {
  tryCatch({
    print(emapplot(pairwise_termsim(gse), showCategory = 50))
  }, error=function(e){warning("Could not plot GO emapplot")})

  tryCatch({
    print(ridgeplot(gse) + labs(x = "enrichment distribution"))
  }, error=function(e){warning("Could not plot GO ridgeplot")})
  write.table(gse@result[,1:11], file=snakemake@output[["gseago_table"]],
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  #gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
}
plot(1, type='n', axes=FALSE, xlab='', ylab='')
dev.off()


## KEGG GSEA
if(nrow(resFilt) > 0){
    print("GSEA on KEGG terms")
    kegg_gene_list <- resFilt$log2FoldChange
    names(kegg_gene_list) <- gene_ids[resFilt$gene]
    kegg_gene_list <- kegg_gene_list[-which(is.na(names(kegg_gene_list)))]
    kegg_gene_list<-sort(na.omit(kegg_gene_list), decreasing = TRUE)

    kk2 <- tryCatch({
        gseKEGG(geneList     = kegg_gene_list,
                 organism     = org_code,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 keyType       = "ncbi-geneid",
                 use_internal_data=TRUE)
        }, error=function(e){NULL})
} else {
    print("No significant genes to run: GSEA on KEGG terms")
    kk2 <- NULL
}

pdf(snakemake@output[["gseakegg_pdf"]], height = 20, width = 20)
if(is.null(gse)){
    write.table(data.frame(), file=snakemake@output[["gseakegg_table"]],
              col.names = TRUE, row.names = TRUE, quote = FALSE)
} else {
    tryCatch({
      print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") +
         facet_grid(.~.sign))
     }, error=function(e){warning("Could not plot KEGG dotplot")})

    tryCatch({
        print(emapplot(pairwise_termsim(kk2), showCategory = 50))
     }, error=function(e){warning("Could not plot KEGG emapplot")})

    tryCatch({
    print(ridgeplot(kk2) + labs(x = "enrichment distribution"))
     }, error=function(e){warning("Could not plot KEGG ridgeplot")})

    write.table(kk2@result[,1:11], file=snakemake@output[["gseakegg_table"]],
              col.names = TRUE, row.names = TRUE, quote = FALSE)
}
plot(1, type='n', axes=FALSE, xlab='', ylab='')
dev.off()
