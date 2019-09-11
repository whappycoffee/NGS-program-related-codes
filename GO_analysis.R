library(BiocGenerics)
library("AnnotationForge")
library(clusterProfiler)
library("org.Hs.eg.db")
library(GOstats)
library(GSEABase)
library(Category)
## iPS stage analysis

iPS = read.table("iPS/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
iPS = na.omit(iPS)
iPS_gene = unique(as.character(iPS$refgene))

gene_iPS = select(org.Hs.eg.db, keys = iPS_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_iPS))
# write.table(na.omit(gene_iPS), file="iPS/iPS_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_iPS$ENTREZID, file="iPS/entrezid_iPS.txt", sep="\t", row.names=FALSE)
#length(gene_iPS$ENTREZID)
iPS_MF = enrichGO(gene = gene_iPS$ENTREZID,
               OrgDb=org.Hs.eg.db,
               ont = "MF",
               pAdjustMethod = "BH",
               minGSSize = 1,
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.01,
               readable = TRUE)
iPS_BP = enrichGO(gene = gene_iPS$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
iPS_CC = enrichGO(gene = gene_iPS$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

iPS_kegg = enrichKEGG(gene = gene_iPS$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(iPS_GO@result)
# write.table(as.data.frame(iPS@result), file = "iPS/CC_iPS.csv")
dotplot(iPS_GO, showCategory =21, title = "MF for iPS")

## S1 stage analysis

S1 = read.table("S1/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S1 = na.omit(S1)
S1_gene = unique(as.character(S1$refgene))

gene_S1 = select(org.Hs.eg.db, keys = S1_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
#dim(na.omit(gene_S1))
# write.table(na.omit(gene_S1), file="S1/S1_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_S1$ENTREZID, file="S1/entrezid_S2.txt", sep="\t", row.names=FALSE)
# length(gene_S1$ENTREZID)
S1_MF = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S1_BP = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S1_CC = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S1_kegg = enrichKEGG(gene = gene_S1$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(S1_GO@result)
# write.table(as.data.frame(S1_GO@result), file = "S1/CC_S1.csv")
dotplot(S1_kegg, showCategory =21, title = "KEGG pathway for S1")

## S2 stage analysis

S2 = read.table("S2/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S2 = na.omit(S2)
S2_gene = unique(as.character(S2$refgene))

gene_S2 = select(org.Hs.eg.db, keys = S2_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_S2))
# write.table(na.omit(gene_S2), file="S2/S2_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_S2$ENTREZID, file="S2/entrezid_S2.txt", sep="\t", row.names=FALSE)
# length(gene_S2$ENTREZID)
S2_MF = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S2_BP = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S2_CC = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S2_kegg = enrichKEGG(gene = gene_S2$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(S2_GO@result)
# write.table(as.data.frame(S2_GO@result), file = "S2/CC_S2.csv")
dotplot(S2_kegg, showCategory =21, title = "KEGG pathway for S2")

## S3 stage analysis

S3 = read.table("S3/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S3 = na.omit(S3)
S3_gene = unique(as.character(S3$refgene))

gene_S3 = select(org.Hs.eg.db, keys = S3_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_S3))
# write.table(na.omit(gene_S3), file="S3/S3_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_iPS$ENTREZID, file="S3/entrezid_S3.txt", sep="\t", row.names=FALSE)
# length(gene_S3$ENTREZID)
S3_MF = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

S3_BP = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S3_CC = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S3_kegg = enrichKEGG(gene = gene_S3$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(S3_GO@result)
# write.table(as.data.frame(S3_GO@result), file = "S3/CC_S3.csv")
dotplot(S3_kegg, showCategory =21, title = "KEGG pathway for S3")

## ISM.D0 stage analysis

ISM.D0 = read.table("ISM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D0 = na.omit(ISM.D0)
ISM.D0_gene = unique(as.character(ISM.D0$refgene))

gene_ISM.D0 = select(org.Hs.eg.db, keys = ISM.D0_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_ISM.D0))
# write.table(na.omit(gene_ISM.D0), file="ISM.D0/ISM.D0_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_ISM.D0$ENTREZID, file="ISM.D0/entrezid_ISM.D0.txt", sep="\t", row.names=FALSE)
# length(gene_ISM.D0$ENTREZID)
ISM.D0_MF = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ISM.D0_BP = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D0_CC = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D0_kegg = enrichKEGG(gene = gene_ISM.D0$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(iPS_GO@result)
# write.table(as.data.frame(ISM.D0_GO@result), file = "ISM.D0/CC_ISM.D0.csv")
dotplot(ISM.D0_kegg, showCategory =21, title = "KEGG pathway for ISM.D0")

## ISM.D4 stage analysis

ISM.D4 = read.table("ISM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D4 = na.omit(ISM.D4)
ISM.D4_gene = unique(as.character(ISM.D4$refgene))

gene_ISM.D4 = select(org.Hs.eg.db, keys = ISM.D4_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_ISM.D4))
# write.table(na.omit(gene_ISM.D4), file="ISM.D4/ISM.D4_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_ISM.D4$ENTREZID, file="ISM.D4/entrezid_ISM.D4.txt", sep="\t", row.names=FALSE)
# length(gene_ISM.D4$ENTREZID)
ISM.D4_MF = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ISM.D4_BP = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D4_CC = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D4_kegg = enrichKEGG(gene = gene_ISM.D4$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(ISM.D4_GO@result)
# write.table(as.data.frame(ISM.D4_GO@result), file = "ISM.D4/CC_ISM.D4.csv")
dotplot(ISM.D4_kegg, showCategory =21, title = "KEGG pathway for ISM.D4")

## ADM.D0 stage analysis

ADM.D0 = read.table("ADM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D0 = na.omit(ADM.D0)
ADM.D0_gene = unique(as.character(ADM.D0$refgene))

gene_ADM.D0 = select(org.Hs.eg.db, keys = ADM.D0_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_ADM.D0))
# write.table(na.omit(gene_ADM.D0), file="ADM.D0/ADM.D0_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_ADM.D0$ENTREZID, file="ADM.D0/entrezid_ADM.D0.txt", sep="\t", row.names=FALSE)
# length(gene_ADM.D0$ENTREZID)
ADM.D0_BP = enrichGO(gene = gene_ADM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ADM.D0_MF = enrichGO(gene = gene_ADM.D0$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
ADM.D0_CC = enrichGO(gene = gene_ADM.D0$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)

ADM.D0_kegg = enrichKEGG(gene = gene_ADM.D0$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(ADM.D0_GO@result)
# write.table(as.data.frame(ADM.D0_GO@result), file = "ADM.D0/CC_ADM.D0.csv")
dotplot(ADM.D0_kegg, showCategory =21, title = "KEGG pathway for ADM.D0")

## ADM.D4 stage analysis

ADM.D4 = read.table("ADM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D4 = na.omit(ADM.D4)
ADM.D4_gene = unique(as.character(ADM.D4$refgene))

gene_ADM.D4 = select(org.Hs.eg.db, keys = ADM.D4_gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
# dim(na.omit(gene_ADM.D4))
# write.table(na.omit(gene_ADM.D4), file="ADM.D4/ADM.D4_entrezid.txt", sep="\t", row.names=FALSE)
# write.table(gene_ADM.D4$ENTREZID, file="ADM.D4/entrezid_ADM.D4.txt", sep="\t", row.names=FALSE)
# length(gene_ADM.D4$ENTREZID)
ADM.D4_BP = enrichGO(gene = gene_ADM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ADM.D4_MF = enrichGO(gene = gene_ADM.D4$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
ADM.D4_CC = enrichGO(gene = gene_ADM.D4$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
ADM.D4_kegg = enrichKEGG(gene = gene_ADM.D4$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
# View(ADM.D4_GO@result)
# write.table(as.data.frame(ADM.D4_GO@result), file = "ADM.D4/CC_ADM.D4.csv")
# View(ADM.D4_kegg@result)
dotplot(ADM.D4_kegg, showCategory =21, title = "KEGG pathway for ADM.D4")
barplot(ADM.D4_GO, showCategory=21,title="Enrichment CC for ADM.D4")
View(ADM.D4_kegg@result)

### creating pie chart
#ggplot(as.data.frame(ADM.D0_GO@result), aes(x = "", y = Count, fill = Description)) + geom_bar(width = 1, stat = "identity")

# counts = ADM.D0_GO@result$Count
# sum(counts)

#### Gostats

enrich <- function(entrezIDs, orgDbName="org.Hs.eg.db", pvalueCutoff=.01){
  require(orgDbName, character.only=TRUE)
  require("GSEABase")
  require("GOstats")
  require("Category")
  goAnn <- get(gsub(".db", "GO", orgDbName))
  universe <- Lkeys(goAnn)
  onto <- c("BP", "MF", "CC")
  res <- lapply(onto, function(.onto){
    param <- new('GOHyperGParams',
                 geneIds= entrezIDs,
                 universeGeneIds=universe,
                 annotation=orgDbName,
                 ontology=.onto,
                 pvalueCutoff=pvalueCutoff,
                 conditional=FALSE,
                 testDirection="over")
    over <- hyperGTest(param)
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(.ids) {
      .sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
      .sym[is.na(.sym)] <- .ids[is.na(.sym)]
      paste(.sym, collapse=";")
    })
    summary <- summary(over)
    if(nrow(summary)>1) summary$Symbols <- glist[as.character(summary[, 1])]
    summary
  })
  names(res) <- onto
  keggAnn <- get(gsub(".db", "PATH", orgDbName))
  universe <- Lkeys(keggAnn)
  param <- new("KEGGHyperGParams",
               geneIds=entrezIDs,
               universeGeneIds=universe,
               annotation=orgDbName,
               categoryName="KEGG",
               pvalueCutoff=pvalueCutoff,
               testDirection="over")
  over <- hyperGTest(param)
  kegg <- summary(over)
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=get(gsub(".db", "SYMBOL", orgDbName)), ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
  })
  kegg$Symbols <- glist[as.character(kegg$KEGGID)]
  res[["kegg"]] <- kegg
  res
}
go = enrich(gene_ADM.D4$ENTREZID,"org.Hs.eg.db", .05)
View(go$MF)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGG.db", version = "3.8")
gene.ADM.D4 = na.omit(gene.ADM.D4)
as.character(gene.ADM.D4)
entrezIDs <- mget(as.character(gene.ADM.D4), org.Hs.egSYMBOL2EG, ifnotfound=NA)
length(entrezIDs)

#### ISM vs ADM vs S2

ISMless3 = read.table("../1-3k_results/ISM.D0/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D0less3 = read.table("../1-3k_results/ADM.D0/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISMless3 = na.omit(ISMless3)
ADM.D0less3 = na.omit(ADM.D0less3)
common_genes1 = unique(intersect(ISMless3$refgene, ADM.D0less3$refgene))
common_genes2 = unique(intersect(ISMless3$refgene, S3$refgene))

gene_ISM_ADM = select(org.Hs.eg.db, keys = common_genes1, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
gene_ISM_S2 = select(org.Hs.eg.db, keys = common_genes2, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ISM_ADM_CC = enrichGO(gene = gene_ISM_ADM$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
dotplot(ISM_ADM_CC, showCategory =21, title = "Cellular Component for ISM.D0 over ADM.D0")+ ggsave("../report/thesis_fig/CC_ISMD0_ADMD0.png")

ISM_ADM_BP = enrichGO(gene = gene_ISM_ADM$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_ADM_BP, showCategory =21, title = "Biological Process for ISM.D0 over ADM.D0")+ ggsave("../report/thesis_fig/BP_ISMD0_ADMD0.png")

ISM_ADM_MF = enrichGO(gene = gene_ISM_ADM$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_ADM_MF, showCategory =21, title = "Molecular Function for ISM.D0 over ADM.D0")+ ggsave("../report/thesis_fig/MF_ISMD0_ADMD0.png")

ISM_ADM_kegg = enrichKEGG(gene = gene_ISM_ADM$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)
dotplot(ISM_ADM_kegg, showCategory =21, title = "KEGG pathway for ISM.D0 over ADM.D0")+ ggsave("../report/thesis_fig/KEGG_ISMD0_ADMD0.png")

ISM_S2_CC = enrichGO(gene = gene_ISM_S2$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_S2_CC, showCategory =21, title = "Cellular Component for ISM.D0 over S2")+ ggsave("../report/thesis_fig/CC_ISMD0_S2.png")

ISM_S2_BP = enrichGO(gene = gene_ISM_S2$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_S2_BP, showCategory =21, title = "Biological Process for ISM.D0 over S2")+ ggsave("../report/thesis_fig/BP_ISMD0_S2.png")

ISM_S2_MF = enrichGO(gene = gene_ISM_S2$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_S2_MF, showCategory =21, title = "Molecular Function for ISM.D0 over S2")+ ggsave("../report/thesis_fig/MF_ISMD0_S2.png")
ISM_S2_kegg = enrichKEGG(gene = gene_ISM_S2$ENTREZID,
                          organism ="human",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          minGSSize = 1,
                          #readable = TRUE ,
                          use_internal_data =FALSE)
dotplot(ISM_S2_kegg, showCategory =21, title = "KEGG pathway for ISM.D0 over S2")+ ggsave("../report/thesis_fig/KEGG_ISMD0_S2.png")

#####  S3 vs ISM vs ADM

ISM_ADM = unique(intersect(ISM.D4$refgene, ADM.D4$refgene))
ISM_S3 = unique(intersect(ISM.D4$refgene, S3$refgene))

gene_ISM_ADM4 = select(org.Hs.eg.db, keys = ISM_ADM, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
gene_ISM_S3 = select(org.Hs.eg.db, keys = ISM_S3, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ISM_ADM4_CC = enrichGO(gene = gene_ISM_ADM4$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_ADM4_CC, showCategory =21, title = "Cellular Component for ISM.D4 over ADM.D4")+ ggsave("../report/thesis_fig/CC_ISMD4_ADMD4.png")
ISM_ADM4_BP = enrichGO(gene = gene_ISM_ADM4$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_ADM4_BP, showCategory =21, title = "Biological Process for ISM.D4 over ADM.D4")+ ggsave("../report/thesis_fig/BP_ISMD4_ADMD4.png")
ISM_ADM4_MF = enrichGO(gene = gene_ISM_ADM4$ENTREZID,
                      OrgDb=org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
dotplot(ISM_ADM4_MF, showCategory =21, title = "Molecular Function for ISM.D4 over ADM.D4")+ ggsave("../report/thesis_fig/MF_ISMD4_ADMD4.png")
ISM_ADM4_kegg = enrichKEGG(gene = gene_ISM_ADM4$ENTREZID,
                          organism ="human",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          minGSSize = 1,
                          #readable = TRUE ,
                          use_internal_data =FALSE)
dotplot(ISM_ADM4_kegg, showCategory =21, title = "KEGG pathway for ISM.D4 over ADM.D4")+ ggsave("../report/thesis_fig/KEGG_ISMD4_ADMD4.png")
ISM_S3_CC = enrichGO(gene = gene_ISM_S3$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
dotplot(ISM_S3_CC, showCategory =21, title = "Cellular Component for ISM.D4 over S3")+ ggsave("../report/thesis_fig/CC_ISMD4_S3.png")
ISM_S3_BP = enrichGO(gene = gene_ISM_S3$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
dotplot(ISM_S3_BP, showCategory =21, title = "Biological Process for ISM.D4 over S3")+ ggsave("../report/thesis_fig/BP_ISMD4_S3.png")
ISM_S3_MF = enrichGO(gene = gene_ISM_S3$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
dotplot(ISM_S3_MF, showCategory =21, title = "Molecular Function for ISM.D4 over S3") + ggsave("../report/thesis_fig/MF_ISMD4_S3.png")
ISM_S3_kegg = enrichKEGG(gene = gene_ISM_S3$ENTREZID,
                         organism ="human",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.01,
                         minGSSize = 1,
                         #readable = TRUE ,
                         use_internal_data =FALSE)

dotplot(ISM_S3_kegg, showCategory =21, title = "KEGG pathway for ISM.D4 over S3") + ggsave("../report/thesis_fig/KEGG_ISMD4_S3.png")

getunique = function(stage,iPS){
  newdata = subset(stage@result[1:21,], !(stage@result[1:21,]$Description %in% iPS@result[1:21,]$Description))
  return(newdata)
}

D4_BP = getunique(ADM.D4_BP,iPS_BP)
D4_MF = getunique(ADM.D4_MF,iPS_MF)
D4_CC = getunique(ADM.D4_CC,iPS_CC)
D4_kegg = getunique(ADM.D4_kegg,iPS_kegg)

### bar chart of go  terms

## iPS line
iPS_go = c(length(iPS_MF@result$ID), length(iPS_BP@result$ID),length(iPS_CC@result$ID),length(iPS_kegg@result$ID))
S1_go = c(length(S1_MF@result$ID), length(S1_BP@result$ID),length(S1_CC@result$ID),length(S1_kegg@result$ID))
S2_go = c(length(S2_MF@result$ID), length(S2_BP@result$ID),length(S2_CC@result$ID),length(S2_kegg@result$ID))
S3_go = c(length(S3_MF@result$ID), length(S3_BP@result$ID),length(S3_CC@result$ID),length(S3_kegg@result$ID))
numbers = c(iPS_go, S1_go, S2_go, S3_go)
terms = rep(c("MF", "BP", "CC", "kegg"),times = 4)
stage = c(rep("iPS", times = 4), rep("S1", times = 4), rep("S2", times = 4), rep("S3", times = 4))
datas = as.data.frame(numbers, terms)
datas$stage = stage
datas$terms = terms
row.names(datas) = seq(1,16)
View(datas)
ggplot(datas, aes(x = terms,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("iPS line GO analysis") + ylab("Numbers") 

a = intersect(iPS_BP@result$Description,S1_BP@result$Description)
b = intersect(S2_BP@result$Description,S3_BP@result$Description)
length(intersect(iPS_BP@result$Description,S1_BP@result$Description))
length(intersect(S2_BP@result$Description,S3_BP@result$Description))
length(intersect(a,b))

a = intersect(iPS_MF@result$Description,S1_MF@result$Description)
b = intersect(S2_MF@result$Description,S3_MF@result$Description)
length(intersect(iPS_MF@result$Description,S1_MF@result$Description))
length(intersect(S2_MF@result$Description,S3_MF@result$Description))
length(intersect(a,b))

a = intersect(iPS_CC@result$Description,S1_CC@result$Description)
b = intersect(S2_CC@result$Description,S3_CC@result$Description)
length(intersect(iPS_CC@result$Description,S1_CC@result$Description))
length(intersect(S2_CC@result$Description,S3_CC@result$Description))
length(intersect(a,b))

a = intersect(iPS_kegg@result$Description,S1_kegg@result$Description)
b = intersect(S2_kegg@result$Description,S3_kegg@result$Description)
length(intersect(iPS_kegg@result$Description,S1_kegg@result$Description))
length(intersect(S2_kegg@result$Description,S3_kegg@result$Description))
length(intersect(a,b))
View(iPS_kegg@result)
### ISM line

ISM.D0_go = c(length(ISM.D0_MF@result$ID), length(ISM.D0_BP@result$ID),length(ISM.D0_CC@result$ID),length(ISM.D0_kegg@result$ID))
ISM.D4_go = c(length(ISM.D4_MF@result$ID), length(ISM.D4_BP@result$ID),length(ISM.D4_CC@result$ID),length(ISM.D4_kegg@result$ID))
numbers = c(ISM.D0_go, ISM.D4_go)
terms = rep(c("MF", "BP", "CC", "kegg"),times = 2)
stage = c(rep("ISM.D0", times = 4), rep("ISM.D4", times = 4))
ISMdatas = as.data.frame(numbers, terms)
ISMdatas$stage = stage
ISMdatas$terms = terms
row.names(ISMdatas) = seq(1,8)
ggplot(ISMdatas, aes(x = terms,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("ISM line GO analysis") + ylab("Numbers") 

### ADM line
ADM.D0_go = c(length(ADM.D0_MF@result$ID), length(ADM.D0_BP@result$ID),length(ADM.D0_CC@result$ID),length(ADM.D0_kegg@result$ID))
ADM.D4_go = c(length(ADM.D4_MF@result$ID), length(ADM.D4_BP@result$ID),length(ADM.D4_CC@result$ID),length(ADM.D4_kegg@result$ID))
numbers = c(ADM.D0_go, ADM.D4_go)
terms = rep(c("MF", "BP", "CC", "kegg"),times = 2)
stage = c(rep("ADM.D0", times = 4), rep("ADM.D4", times = 4))
ADMdatas = as.data.frame(numbers, terms)
ADMdatas$stage = stage
ADMdatas$terms = terms
row.names(ADMdatas) = seq(1,8)
ggplot(ADMdatas, aes(x = terms,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("ADM line GO analysis") + ylab("Numbers") 


####### functional analysis for all shared genes

gene_all = select(org.Hs.eg.db, keys = g, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
all_MF = enrichGO(gene = gene_all$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
all_BP = enrichGO(gene = gene_all$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
all_CC = enrichGO(gene = gene_all$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

all_kegg = enrichKEGG(gene = gene_all$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

numbers = c(length(all_MF@result$ID), length(all_BP@result$ID),length(all_CC@result$ID),length(all_kegg@result$ID))
all = as.data.frame(numbers)
all$terms = c("MF", "BP", "CC", "kegg")
View(all)
ggplot(all, aes(x = terms,y = numbers)) + geom_bar(stat = "identity") + 
  xlab("all shared genes GO analysis") + ylab("Numbers")

###### funcitonal analysis for active genes

gene_iPS = select(org.Hs.eg.db, keys = as.character(gene.iPS), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
iPS_MF = enrichGO(gene = gene_iPS$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
iPS_BP = enrichGO(gene = gene_iPS$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
iPS_CC = enrichGO(gene = gene_iPS$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

iPS_kegg = enrichKEGG(gene = gene_iPS$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_S1 = select(org.Hs.eg.db, keys = as.character(gene.S1), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
S1_MF = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S1_BP = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S1_CC = enrichGO(gene = gene_S1$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

S1_kegg = enrichKEGG(gene = gene_S1$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_S2 = select(org.Hs.eg.db, keys = as.character(gene.S2), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
S2_MF = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S2_BP = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S2_CC = enrichGO(gene = gene_S2$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

S2_kegg = enrichKEGG(gene = gene_S2$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_S3 = select(org.Hs.eg.db, keys = as.character(gene.S3), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
S3_MF = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S3_BP = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
S3_CC = enrichGO(gene = gene_S3$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

S3_kegg = enrichKEGG(gene = gene_S3$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_ISM.D0 = select(org.Hs.eg.db, keys = as.character(gene.ISM.D0), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ISM.D0_MF = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D0_BP = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D0_CC = enrichGO(gene = gene_ISM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ISM.D0_kegg = enrichKEGG(gene = gene_ISM.D0$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_ISM.D4 = select(org.Hs.eg.db, keys = as.character(gene.ISM.D4), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ISM.D4_MF = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D4_BP = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ISM.D4_CC = enrichGO(gene = gene_ISM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ISM.D4_kegg = enrichKEGG(gene = gene_ISM.D4$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_ADM.D0 = select(org.Hs.eg.db, keys = as.character(gene.ADM.D0), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ADM.D0_MF = enrichGO(gene = gene_ADM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ADM.D0_BP = enrichGO(gene = gene_ADM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ADM.D0_CC = enrichGO(gene = gene_ADM.D0$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ADM.D0_kegg = enrichKEGG(gene = gene_ADM.D0$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)

gene_ADM.D4 = select(org.Hs.eg.db, keys = as.character(gene.ADM.D4), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ADM.D4_MF = enrichGO(gene = gene_ADM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ADM.D4_BP = enrichGO(gene = gene_ADM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)
ADM.D4_CC = enrichGO(gene = gene_ADM.D4$ENTREZID,
                  OrgDb=org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01,
                  readable = TRUE)

ADM.D4_kegg = enrichKEGG(gene = gene_ADM.D4$ENTREZID,
                      organism ="human",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      minGSSize = 1,
                      #readable = TRUE ,
                      use_internal_data =FALSE)



iPS_go = c(length(iPS_MF@result$ID), length(iPS_BP@result$ID),length(iPS_CC@result$ID),length(iPS_kegg@result$ID))
S1_go = c(length(S1_MF@result$ID), length(S1_BP@result$ID),length(S1_CC@result$ID),length(S1_kegg@result$ID))
S2_go = c(length(S2_MF@result$ID), length(S2_BP@result$ID),length(S2_CC@result$ID),length(S2_kegg@result$ID))
S3_go = c(length(S3_MF@result$ID), length(S3_BP@result$ID),length(S3_CC@result$ID),length(S3_kegg@result$ID))
numbers = c(iPS_go, S1_go, S2_go, S3_go)
terms = rep(c("MF", "BP", "CC", "kegg"),times = 4)
stage = c(rep("iPS", times = 4), rep("S1", times = 4), rep("S2", times = 4), rep("S3", times = 4))
datas = as.data.frame(numbers, terms)
datas$stage = stage
datas$terms = terms
row.names(datas) = seq(1,16)
View(datas)
ggplot(datas, aes(x = terms,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("iPS line GO analysis") + ylab("Numbers") 

a = intersect(ADM.D0_BP@result$Description,ADM.D4_BP@result$Description)
length(a)
b = intersect(S2_BP@result$Description,S3_BP@result$Description)
length(intersect(iPS_BP@result$Description,S1_BP@result$Description))
length(intersect(S2_BP@result$Description,S3_BP@result$Description))
length(intersect(a,b))

a = intersect(iPS_MF@result$Description,S1_MF@result$Description)
b = intersect(S2_MF@result$Description,S3_MF@result$Description)
length(intersect(iPS_MF@result$Description,S1_MF@result$Description))
length(intersect(S2_MF@result$Description,S3_MF@result$Description))
length(intersect(a,b))

a = intersect(iPS_CC@result$Description,S1_CC@result$Description)
b = intersect(S2_CC@result$Description,S3_CC@result$Description)
length(intersect(iPS_CC@result$Description,S1_CC@result$Description))
length(intersect(S2_CC@result$Description,S3_CC@result$Description))
length(intersect(a,b))

a = intersect(iPS_kegg@result$Description,S1_kegg@result$Description)
b = intersect(S2_kegg@result$Description,S3_kegg@result$Description)
length(intersect(iPS_kegg@result$Description,S1_kegg@result$Description))
length(intersect(S2_kegg@result$Description,S3_kegg@result$Description))
length(intersect(a,b))
dotplot(iPS_MF, showCategory =21, title = "MF for iPS")
dotplot(iPS_BP, showCategory =21, title = "BP for iPS")
dotplot(iPS_CC, showCategory =21, title = "CC for iPS")
dotplot(iPS_kegg, showCategory =21, title = "kegg for iPS")
dotplot(S1_MF, showCategory =21, title = "MF for S1")
dotplot(S1_BP, showCategory =21, title = "BP for S1")
dotplot(S1_CC, showCategory =21, title = "CC for S1")
dotplot(S1_kegg, showCategory =21, title = "kegg for S1")
dotplot(S2_MF, showCategory =21, title = "MF for S2")
dotplot(S2_BP, showCategory =21, title = "BP for S2")
dotplot(S2_CC, showCategory =21, title = "CC for S2")
dotplot(S2_kegg, showCategory =21, title = "kegg for S2")
dotplot(S3_MF, showCategory =21, title = "MF for S3")
dotplot(S3_BP, showCategory =21, title = "BP for S3")
dotplot(S3_CC, showCategory =21, title = "CC for S3")
dotplot(S3_kegg, showCategory =21, title = "kegg for S3")
dotplot(ISM.D0_MF, showCategory =21, title = "MF for ISM.D0")
dotplot(ISM.D0_BP, showCategory =21, title = "BP for ISM.D0")
dotplot(ISM.D0_CC, showCategory =21, title = "CC for ISM.D0")
dotplot(ISM.D0_kegg, showCategory =21, title = "kegg for ISM.D0")
dotplot(ISM.D4_MF, showCategory =21, title = "MF for ISM.D4")
dotplot(ISM.D4_BP, showCategory =21, title = "BP for ISM.D4")
dotplot(ISM.D4_CC, showCategory =21, title = "CC for ISM.D4")
dotplot(ISM.D4_kegg, showCategory =21, title = "kegg for ISM.D4")
dotplot(ADM.D0_MF, showCategory =21, title = "MF for ADM.D0")
dotplot(ADM.D0_BP, showCategory =21, title = "BP for ADM.D0")
dotplot(ADM.D0_CC, showCategory =21, title = "CC for ADM.D0")
dotplot(ADM.D0_kegg, showCategory =21, title = "kegg for ADM.D0")
dotplot(ADM.D4_MF, showCategory =21, title = "MF for ADM.D4")
dotplot(ADM.D4_BP, showCategory =21, title = "BP for ADM.D4")
dotplot(ADM.D4_CC, showCategory =21, title = "CC for ADM.D4")
dotplot(ADM.D4_kegg, showCategory =21, title = "kegg for ADM.D4")



