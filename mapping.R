library(refGenome)
iso = ensemblGenome()
ips = read.csv("Documents/Research/tapis_run/iPS_tapis_out/AS_summary.csv", header = TRUE,sep = ",")
geneid = ips$Gene
read.gtf(iso, filename="Documents/Research/tapis_run/iPS_tapis_out/assembled.gtf")
genes = gtf@ev$gtf[ ,c("gene_name","gene_id","transcript_id","feature", "")]

genenames = gtf@ev[["gtf"]][["gene_name"]]
transname = gtf@ev[["gtf"]][["transcript_name"]]
chr = iso@ev[["gtf"]][["seqid"]]
chr = chr[1:50]
star  = iso@ev[["gtf"]][["start"]]
star = star[1:50]
end = iso@ev[["gtf"]][["end"]]
end = end[1:50]
strand = iso@ev[["gtf"]][["strand"]]
strand = strand[1:50]
table(genenames)
length(table(genenames))
length(table(transname))
table(transID)
transID = gtf@ev[["gtf"]][["transcript_id"]]
table(transID)
length(table(transID))
hg19 = ensemblGenome()
read.gtf(hg19, filename="Documents/Research/data/Homo_sapiens.GRCh37.87.gtf")

library(biomaRt)
library(Gviz)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X") 
options(ucscChromosomeNames = FALSE)
gnm <- GRanges("1:23409672-23409850")
gat <- GenomeAxisTrack(range = gnm)
gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))
gtx <- GeneRegionTrack(gnm_gns, name = "tx", geneSymbol = TRUE,
                       showId = TRUE)
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm)
plotTracks(list(ht))
genomeToTranscript(gnm, edbx)

listMarts()
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
filters = listFilters(grch37)[1:3,]
filters[1:3,]
attributes = listAttributes(ensembl)
listAttributes(grch37)
attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "ensembl_peptide_id", "chromosome_name","start_position", "transcript_start","transcript_end","end_position","ccds","rank", "hgnc_symbol")
getBM(attributes = attributes, filters = c("chromosome_name", "start","end"), value = list(chr, star, end),mart = grch37)

region = data.frame(chr,start,end, strand)
length(start)
star
result = getBM(attributes = attributes, filters = c("chromosome_name", "start","end"), value = list(chr[1:2], star[1:2], end[1:2]),mart = grch37)
result
View(result)

getid = getBM(attributes = attributes, filters = "ensembl_gene_id", value = geneid,mart = grch37)
View(getid)
