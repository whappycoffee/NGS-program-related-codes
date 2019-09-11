library(edgeR)
group = c("iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D4", "ADM.D0")
libSizes <- as.vector(colSums(data))
d <- DGEList(counts=data,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.4
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="outputFile.txt",sep="\t")
View(as.matrix(results$table))
### iPS ###
iPS_S1 = data[c("iPS", "S1")]
libSizes <- as.vector(colSums(iPS_S1))
group = c("iPS", "S1")
d <- DGEList(counts=iPS_S1,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_S1.txt",sep="\t")
View(as.matrix(results$table))

iPS_S2 = data[c("iPS", "S2")]
libSizes <- as.vector(colSums(iPS_S2))
group = c("iPS", "S2")
d <- DGEList(counts=iPS_S2,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_S2.txt",sep="\t")
View(as.matrix(results$table))

iPS_S3 = data[c("iPS", "S3")]
libSizes <- as.vector(colSums(iPS_S3))
group = c("iPS", "S3")
d <- DGEList(counts=iPS_S3,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_S3.txt",sep="\t")
View(as.matrix(results$table))

iPS_ISM.D0 = data[c("iPS", "ISM.D0")]
libSizes <- as.vector(colSums(iPS_ISM.D0))
group = c("iPS", "ISM.D0")
d <- DGEList(counts=iPS_ISM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_ISM.D0.txt",sep="\t")
View(as.matrix(results$table))

iPS_ISM.D4 = data[c("iPS", "ISM.D4")]
libSizes <- as.vector(colSums(iPS_ISM.D4))
group = c("iPS", "ISM.D4")
d <- DGEList(counts=iPS_ISM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_ISM.D4.txt",sep="\t")
View(as.matrix(results$table))

iPS_ADM.D0 = data[c("iPS", "ADM.D0")]
libSizes <- as.vector(colSums(iPS_ADM.D0))
group = c("iPS", "ADM.D0")
d <- DGEList(counts=iPS_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_ADM.D0.txt",sep="\t")
View(as.matrix(results$table))

iPS_ADM.D4 = data[c("iPS", "ADM.D4")]
libSizes <- as.vector(colSums(iPS_ADM.D4))
group = c("iPS", "ADM.D4")
d <- DGEList(counts=iPS_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="iPS_ADM.D4.txt",sep="\t")
View(as.matrix(results$table))

## S1

S1_S2 = data[c("S1", "S2")]
libSizes <- as.vector(colSums(S1_S2))
group = c("S1", "S2")
d <- DGEList(counts=S1_S2,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_S2.txt",sep="\t")
View(as.matrix(results$table))

S1_S3 = data[c("S1", "S3")]
libSizes <- as.vector(colSums(S1_S3))
group = c("S1", "S3")
d <- DGEList(counts=S1_S3,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_S3.txt",sep="\t")
View(as.matrix(results$table))

S1_ISM.D0 = data[c("S1", "ISM.D0")]
libSizes <- as.vector(colSums(S1_ISM.D0))
group = c("S1", "ISM.D0")
d <- DGEList(counts=S1_ISM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_ISM.D0.txt",sep="\t")
View(as.matrix(results$table))

S1_ISM.D4 = data[c("S1", "ISM.D4")]
libSizes <- as.vector(colSums(S1_ISM.D4))
group = c("S1", "ISM.D4")
d <- DGEList(counts=S1_ISM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_ISM.D4.txt",sep="\t")
View(as.matrix(results$table))

S1_ADM.D0 = data[c("S1", "ADM.D0")]
libSizes <- as.vector(colSums(S1_ADM.D0))
group = c("S1", "ADM.D0")
d <- DGEList(counts=S1_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_ADM.D0.txt",sep="\t")
View(as.matrix(results$table))

S1_ADM.D4 = data[c("S1", "ADM.D4")]
libSizes <- as.vector(colSums(S1_ADM.D4))
group = c("S1", "ADM.D4")
d <- DGEList(counts=S1_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S1_ADM.D4.txt",sep="\t")
View(as.matrix(results$table))

#### S2

S2_S3 = data[c("S2", "S3")]
libSizes <- as.vector(colSums(S2_S3))
group = c("S2", "S3")
d <- DGEList(counts=S2_S3,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S2_S3.txt",sep="\t")
View(as.matrix(results$table))

S2_ISM.D0 = data[c("S2", "ISM.D0")]
libSizes <- as.vector(colSums(S2_ISM.D0))
group = c("S2", "ISM.D0")
d <- DGEList(counts=S2_ISM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S2_ISM.D0.txt",sep="\t")
View(as.matrix(results$table))

S2_ISM.D4 = data[c("S2", "ISM.D4")]
libSizes <- as.vector(colSums(S2_ISM.D4))
group = c("S2", "ISM.D4")
d <- DGEList(counts=S2_ISM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S2_ISM.D4.txt",sep="\t")
View(as.matrix(results$table))

S2_ADM.D0 = data[c("S2", "ADM.D0")]
libSizes <- as.vector(colSums(S2_ADM.D0))
group = c("S2", "ADM.D0")
d <- DGEList(counts=S2_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S2_ADM.D0.txt",sep="\t")
View(as.matrix(results$table))

S2_ADM.D4 = data[c("S2", "ADM.D4")]
libSizes <- as.vector(colSums(S2_ADM.D4))
group = c("S2", "ADM.D4")
d <- DGEList(counts=S2_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S2_ADM.D4.txt",sep="\t")
View(as.matrix(results$table))


#### S3


S3_ISM.D0 = data[c("S3", "ISM.D0")]
libSizes <- as.vector(colSums(S3_ISM.D0))
group = c("S3", "ISM.D0")
d <- DGEList(counts=S3_ISM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S3_ISM.D0.txt",sep="\t")
View(as.matrix(results$table))

S3_ISM.D4 = data[c("S3", "ISM.D4")]
libSizes <- as.vector(colSums(S3_ISM.D4))
group = c("S3", "ISM.D4")
d <- DGEList(counts=S3_ISM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S3_ISM.D4.txt",sep="\t")

S3_ADM.D0 = data[c("S3", "ADM.D0")]
libSizes <- as.vector(colSums(S3_ADM.D0))
group = c("S3", "ADM.D0")
d <- DGEList(counts=S3_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S3_ADM.D0.txt",sep="\t")

S3_ADM.D4 = data[c("S3", "ADM.D4")]
libSizes <- as.vector(colSums(S3_ADM.D4))
group = c("S3", "ADM.D4")
d <- DGEList(counts=S3_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="S3_ADM.D4.txt",sep="\t")

#### ISM.D0

ISM.D0_ISM.D4 = data[c("ISM.D0", "ISM.D4")]
libSizes <- as.vector(colSums(ISM.D0_ISM.D4))
group = c("ISM.D0", "ISM.D4")
d <- DGEList(counts=ISM.D0_ISM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ISM.D0_ISM.D4.txt",sep="\t")

ISM.D0_ADM.D0 = data[c("ISM.D0", "ADM.D0")]
libSizes <- as.vector(colSums(ISM.D0_ADM.D0))
group = c("ISM.D0", "ADM.D0")
d <- DGEList(counts=ISM.D0_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ISM.D0_ADM.D0.txt",sep="\t")

ISM.D0_ADM.D4 = data[c("ISM.D0", "ADM.D4")]
libSizes <- as.vector(colSums(ISM.D0_ADM.D4))
group = c("ISM.D0", "ADM.D4")
d <- DGEList(counts=ISM.D0_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ISM.D0_ADM.D4.txt",sep="\t")

### ISM.D4

ISM.D4_ADM.D0 = data[c("ISM.D4", "ADM.D0")]
libSizes <- as.vector(colSums(ISM.D4_ADM.D0))
group = c("ISM.D4", "ADM.D0")
d <- DGEList(counts=ISM.D4_ADM.D0,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ISM.D4_ADM.D0.txt",sep="\t")

ISM.D4_ADM.D4 = data[c("ISM.D4", "ADM.D4")]
libSizes <- as.vector(colSums(ISM.D4_ADM.D4))
group = c("ISM.D4", "ADM.D4")
d <- DGEList(counts=ISM.D4_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ISM.D4_ADM.D4.txt",sep="\t")

#### ADM.D0

ADM.D0_ADM.D4 = data[c("ADM.D0", "ADM.D4")]
libSizes <- as.vector(colSums(ADM.D0_ADM.D4))
group = c("ADM.D0", "ADM.D4")
d <- DGEList(counts=ADM.D0_ADM.D4,group=group,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 0.05
de.com <- exactTest(d)
results <- topTags(de.com,n = length(data[,1]))
write.table(as.matrix(results$table),file="ADM.D0_ADM.D4.txt",sep="\t")
