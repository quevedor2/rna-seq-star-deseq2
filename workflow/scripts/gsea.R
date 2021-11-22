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
library("msigdbr")


# Create a reference map of ENSEMBL to SYMBOL
if(snakemake@params[['genome']] == 'mus_musculus'){
    library("org.Mm.eg.db")
    annotation <- "org.Mm.eg.db"
    species <- 'Mus musculus'
    genome <- org.Mm.eg.db
    org_code <- 'mmu'
} else if(snakemake@params[['genome']] == 'homo_sapiens'){
    library("org.Hs.eg.db")
    annotation <- "org.Hs.eg.db"
    species <- 'Homo sapiens'
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
resSig = res[ which(res$padj < as.numeric(as.character(snakemake@params[["maxp"]])) &
                         abs(res$log2FoldChange) > 0 &
                         res$baseMean > as.numeric(as.character(snakemake@params[["minbase"]]))), ]
resFilt <- res[res$baseMean > as.numeric(as.character(snakemake@params[["minbase"]])),]


msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C1'=list(NULL),                      # positional gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                  'C6'=list(NULL),                      # oncogenic signature gene sets
                  'C7'=list('IMMUNESIGDB'),             # immunologic signature gene sets
                  'C8'=list(NULL))                      # cell type signature gene sets
ora_gra <- lapply(names(msig_lvls), function(mlvl){
  sub_ora_gra <- lapply(msig_lvls[[mlvl]], function(sublvl){
    print(paste0(">", mlvl, ":", sublvl, "..."))
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
      dplyr::select(gs_name, entrez_gene) %>% 
      as.data.frame()
    
    # overrepresentation analysis
    sig_ora <- tryCatch({
      enricher(gene = na.omit(gene_ids[resSig$gene]), TERM2GENE = msig_ds)@result
    }, error=function(e){NULL})
    
    
    
    # GSEA analysis
    lfc_v <- setNames(resFilt$log2FoldChange, 
                      gene_ids[resFilt$gene])
    msig_gsea <- tryCatch({
      GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, pvalueCutoff = 1)
    })
    msig_gsea_df <- as.data.frame(msig_gsea)[,1:10]
    
    
    ## network of similar terms
    n <- sum(msig_gsea@result$pvalue < 0.05)
    if(n > 25) n <- 25
    if(n < 5) n <- 5
    emap <- emapplot(pairwise_termsim(msig_gsea), showCategory = n)
    
    ## Dotplot of upregulate and downregulated terms
    dp <- dotplot(msig_gsea, showCategory = n, title = "Enriched Pathways" , split=".sign") +
      facet_grid(.~.sign)
    
    ## GSEA-type plot of terms
    gg_gsea <- gseaplot2(msig_gsea, geneSetID = head(msig_gsea@result$ID,n))
    
    ## GSEA-type plot of terms
    ridge <- ridgeplot(msig_gsea, showCategory = n) + labs(x = "enrichment distribution")
    
    
    return(list("gsea-viz"=list(emap, dp, gg_gsea, ridge),
                "gsea-tbl"=msig_gsea_df,
                "ora-tbl"=sig_ora))
  })
  return(sub_ora_gra)
})
ora_gra <- unlist(ora_gra, recursive = F)
names(ora_gra) <- unlist(sapply(names(msig_lvls), function(id) paste(id, msig_lvls[[id]], sep="-")))

# Plot out all levels of msigdb
pdf(snakemake@output[["gsea_pdf"]], height = 15, width = 15)
lapply(names(ora_gra), function(mlvl){
  sapply(ora_gra[[mlvl]][['gsea-viz']], function(igv) {
    try(print(igv + ggtitle(mlvl)), silent=T)
  })
  plot(1, type='n', axes=FALSE, xlab='', ylab='')
})
dev.off()

# Write out the GSEA table
gsea_tbl <- do.call(rbind, lapply(ora_gra, function(i) i$`gsea-tbl`))
gsea_tbl <- as.data.frame(gsea_tbl) 
gsea_tbl$msig_lvl <- gsub("^(.*?)\\..*", "\\1", rownames(gsea_tbl))
write.table(gsea_tbl, file=snakemake@output[["gsea_tbl"]],
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Write out the over-representation analysis table
ora_tbl <- do.call(rbind, lapply(ora_gra, function(i) i$`ora-tbl`))
ora_tbl <- as.data.frame(ora_tbl) 
ora_tbl$msig_lvl <- gsub("^(.*?)\\..*", "\\1", rownames(ora_tbl))
write.table(ora_tbl, file=snakemake@output[["ora_tbl"]],
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
