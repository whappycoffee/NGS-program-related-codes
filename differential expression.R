library("edgeR")
library(plyr)
data_raw <- read.table("GSE49712_HTSeq.txt.gz", header = TRUE)
iPS = read.table("iPS/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S1 = read.table("S1/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S2 = read.table("S2/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S3 = read.table("S3/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D0 = read.table("ISM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D4 = read.table("ISM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D0 = read.table("ADM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D4 = read.table("ADM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )

iPS = as.data.frame(table(na.omit(iPS$refgene)))
colnames(iPS) = c("gene", "iPS")
S1 = as.data.frame(table(na.omit(S1$refgene)))
colnames(S1) = c("gene", "S1")
S2 = as.data.frame(table(na.omit(S2$refgene)))
colnames(S2) = c("gene", "S2")
S3 = as.data.frame(table(na.omit(S3$refgene)))
colnames(S3) = c("gene", "S3")
ISM.D0 = as.data.frame(table(na.omit(ISM.D0$refgene)))
colnames(ISM.D0) = c("gene", "ISM.D0")
ISM.D4 = as.data.frame(table(na.omit(ISM.D4$refgene)))
colnames(ISM.D4) = c("gene", "ISM.D4")
ADM.D0 = as.data.frame(table(na.omit(ADM.D0$refgene)))
colnames(ADM.D0) = c("gene", "ADM.D0")
ADM.D4 = as.data.frame(table(na.omit(ADM.D4$refgene)))
colnames(ADM.D4) = c("gene", "ADM.D4")

data = Reduce(function(x, y) merge(x, y, all=TRUE), list(iPS, S1, S2,S3,ISM.D0,ISM.D4,ADM.D0,ADM.D4))
data[is.na(data)] = 0
row.names(data) = data$gene
data[,1] = NULL
dim(data)

cpm_log <- cpm(data, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data <- data[median_log2_cpm > expr_cutoff, ]
dim(data)
dim(cpm_log)
head(cpm_log)
heatmap(cor(cpm_log))
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
summary(pca)
group = c("iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4")
group1 = c("iPS", "progenitor", "myoblasts", "myotubes", "myoblasts", "myotubes", "myoblasts", "myotubes")
group2 = c("iPS", "S1", "S", "S", "ISM", "ISM", "ADM", "ADM")
y1 <- DGEList(counts = data, group = group1)
y2 <- DGEList(counts = data, group = group2)
y <- DGEList(counts = data, group = group)
y
y <- calcNormFactors(y, method = "TMM")
y1 <- calcNormFactors(y1, method = "TMM")
y2 <- calcNormFactors(y2, method = "TMM")

design <- model.matrix(~group)
design1 <- model.matrix(~group1)
design2 <- model.matrix(~group2)

y <- estimateDisp(y)
sqrt(y$common.dispersion)
y1 <- estimateDisp(y1)
y2 <- estimateDisp(y2)
sqrt(y1$common.dispersion)
sqrt(y2$common.dispersion)

plotBCV(y1)
plotBCV(y2)

et1 <- exactTest(y1)
results_edgeR1 <- topTags(et1, n = nrow(data), sort.by = "none")
head(results_edgeR1$table)

et2 <- exactTest(y2)
results_edgeR2 <- topTags(et2, n = nrow(data), sort.by = "none")
head(results_edgeR2$table)

sum(results_edgeR1$table$FDR < .1)
sum(results_edgeR2$table$FDR < .1)

plotSmear(et1, de.tags = rownames(results_edgeR1)[results_edgeR1$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")

boxplot(as.numeric(data["DMD", ]) ~ group2)

y <- estimateDisp(y,design)

fit <- glmFit(y, design)
y
plotMDS(y)
y$samples
plotBCV(y)

#### ADM vs ISM

ADM_ISM = data[c("ISM.D0", "ADM.D0")]
group = factor(c("ISM", "ADM"))
y_ADM_ISM <- DGEList(counts=ADM_ISM,group=group)
y_ADM_ISM <- calcNormFactors(y_ADM_ISM)
design <- model.matrix(~group)
y_ADM_ISM <- estimateDisp(y_ADM_ISM,design)
plotBCV(y_ADM_ISM)
fit <- glmFit(y_ADM_ISM, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

et_ADM_ISM <- exactTest(y_ADM_ISM)
results_edgeR_ADM_ISM <- topTags(et_ADM_ISM, n = nrow(ADM_ISM), sort.by = "none")

sum(results_edgeR_ADM_ISM$table$FDR < .1)

plotSmear(et_ADM_ISM, de.tags = rownames(results_edgeR_ADM_ISM)[results_edgeR_ADM_ISM$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")

### MDS plot
library(magrittr)
library(dplyr)
library(ggpubr)

mds <- data %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = rownames(swiss),
          size = 1,
          repel = TRUE)

library(Rsubread)

counts =featureCounts(files="iPS/ToFU/hq_isoforms.fastq.sorted.sam",annot.ext="../data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf",
                      isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")
View(counts$stat)
table = processAmplicons(readfile = "iPS/ToFU/hq_isoforms.fastq", )

### no replica
keep <- rowSums(cpm(y)>1) >= 2
y_keep = y[keep, , keep.lib.sizes=FALSE]
dim(y_keep$samples)
View(y_keep$samples)
View(y_keep$counts)
y_norm = calcNormFactors(y_keep)

bcv <- 0.4
counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2)
et <- exactTest(y_norm, dispersion=bcv^2)
y_norm$samples$group <- 1
y_norm$samples$group
y0 <- estimateDisp(y_norm[housekeeping,], trend="none", tagwise=FALSE)
