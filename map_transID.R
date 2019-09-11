library(refGenome)
library(plyr)
library(biomaRt)
library(Gviz)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
ADM.D0 = ensemblGenome()
ADM.D4 = ensemblGenome()
ISM.D0 = ensemblGenome()
ISM.D4 = ensemblGenome()
S1 = ensemblGenome()
S2 = ensemblGenome()
S3 = ensemblGenome()
iPS = ensemblGenome()
## read GTF
read.gtf(iPS,"iPS/tapis_out/assembled.gtf")
read.gtf(ADM.D0, "ADM.D0/tapis_out/assembled.gtf")
read.gtf(ADM.D4, "ADM.D4/tapis_out/assembled.gtf")
read.gtf(ISM.D0, "ISM.D0/tapis_out/assembled.gtf")
read.gtf(ISM.D4, "ISM.D4/tapis_out/assembled.gtf")
read.gtf(S1, "S1/tapis_out/assembled.gtf")
read.gtf(S2, "S2/tapis_out/assembled.gtf")
read.gtf(S3, "S3/tapis_out/assembled.gtf")

gettrnasid = function(stage){
  info = stage@ev$gtf[c("seqid", "start", "end", "transcript_id")]
  newdata = ddply(info, .(transcript_id, seqid), function(x) c( min(x[,2]),max(x[,3])))
}

map = gettrnasid(ISM.D4)
View(map)
chr = map$seqid
star = map$V1
end = map$V2
## maping with ensemble

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
filters = listFilters(grch37)[1:3,]
listAttributes(grch37)
attributes =c("ensembl_transcript_id","transcript_start", "transcript_end")
mapresult = getBM(attributes = attributes, filters = c("chromosome_name", "start","end"), value = list(chr, star, end),mart = grch37)
View(mapresult)
