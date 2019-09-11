ADM.D0 = read.table("iPS/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S1 = read.table("S1/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S2 = read.table("S2/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
S3 = read.table("S3/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D0 = read.table("ISM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ISM.D4 = read.table("ISM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D0 = read.table("ADM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )
ADM.D4 = read.table("ADM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t", header = T )

## iPS stage
matched = na.omit(iPS)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(iPS$pbgene))
iPS_gene_not_matched = total_gene - gene_matched
iPS_iso_not_matched = length(iPS$pbid) - length(nona_iPS$pbid)

## S1 stage
matched = na.omit(S1)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(S1$pbgene))
S1_gene_not_matched = total_gene - gene_matched
S1_iso_not_matched = length(S1$pbid) - length(matched$pbid)

## S2 stage
matched = na.omit(S2)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(S2$pbgene))
S2_gene_not_matched = total_gene - gene_matched
S2_iso_not_matched = length(S2$pbid) - length(matched$pbid)

## S2 stage
matched = na.omit(S2)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(S2$pbgene))
S2_gene_not_matched = total_gene - gene_matched
S2_iso_not_matched = length(S2$pbid) - length(matched$pbid)

## S3 stage
matched = na.omit(S3)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(S3$pbgene))
S3_gene_not_matched = total_gene - gene_matched
S3_iso_not_matched = length(S3$pbid) - length(matched$pbid)

## ISM.D0 stage
matched = na.omit(ISM.D0)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(ISM.D0$pbgene))
ISM.D0_gene_not_matched = total_gene - gene_matched
ISM.D0_iso_not_matched = length(ISM.D0$pbid) - length(matched$pbid)

## ISM.D4 stage
matched = na.omit(ISM.D4)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(ISM.D4$pbgene))
ISM.D4_gene_not_matched = total_gene - gene_matched
ISM.D4_iso_not_matched = length(ISM.D4$pbid) - length(matched$pbid)

## ADM.D0 stage
matched = na.omit(ADM.D0)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(ADM.D0$pbgene))
ADM.D0_gene_not_matched = total_gene - gene_matched
ADM.D0_iso_not_matched = length(ADM.D0$pbid) - length(matched$pbid)

## ADM.D4 stage
matched = na.omit(ADM.D4)
gene_matched = length(unique(matched$pbgene))
total_gene = length(unique(ADM.D4$pbgene))
ADM.D4_gene_not_matched = total_gene - gene_matched
ADM.D4_iso_not_matched = length(ADM.D4$pbid) - length(matched$pbid)

## creating data frame

numbers = c(iPS_gene_not_matched, iPS_iso_not_matched, 
            S1_gene_not_matched, S1_iso_not_matched,
            S2_gene_not_matched, S2_iso_not_matched,
            S3_gene_not_matched, S3_iso_not_matched,
            ISM.D0_gene_not_matched, ISM.D0_iso_not_matched,
            ISM.D4_gene_not_matched, ISM.D4_iso_not_matched,
            ADM.D0_gene_not_matched, ADM.D0_iso_not_matched,
            ADM.D4_gene_not_matched, ADM.D4_iso_not_matched)
type = rep(c("gene", "isoform"), times = 8)

stage = c("iPS", "iPS",
          "S1", "S1",
          "S2", "S2",
          "S3", "S3",
          "ISM.D0", "ISM.D0",
          "ISM.D4", "ISM.D4",
          "ADM.D0", "ADM.D0",
          "ADM.D4", "ADM.D4")
novel = as.data.frame(numbers, type)
novel$type = type
novel$stage = stage
row.names(novel) = seq(1,16)
View(novel)
Stage = factor(stage, levels = c("iPS","S1","S2","S3", "ISM.D0", "ISM.D4", "ADM.D0","ADM.D4"))

ggplot(novel, aes(x = Stage, y = numbers, fill = type)) + geom_bar(stat = "identity", position = "dodge") + xlab("stages") + ylab("numbers")


