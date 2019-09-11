library(ggplot2)
require(venneuler)
library(VennDiagram)
iPS = read.delim("~/Documents/Research/results/iPS/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
ADM.D0 = read.delim("~/Documents/Research/results/ADM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
ADM.D4 = read.delim("~/Documents/Research/results/ADM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
ISM.D0 = read.delim("~/Documents/Research/results/ISM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
ISM.D4 = read.delim("~/Documents/Research/results/ISM.D4/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
S1 = read.delim("~/Documents/Research/results/S1/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
S2 = read.delim("~/Documents/Research/results/S2/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
S3 = read.delim("~/Documents/Research/results/S3/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")

ADM.D0less3 = read.delim("~/Documents/Research/1-3k_results/ADM.D0/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
ISM.D0less3 = read.delim("~/Documents/Research/1-3k_results/ISM.D0//iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")

ADM.D0less3 = na.omit(ADM.D0less3)
ISM.D0less3 = na.omit(ISM.D0less3)
iPS = na.omit(iPS)
S1 = na.omit(S1)
S2 = na.omit(S2)
S3 = na.omit(S3)
ISM.D0 = na.omit(ISM.D0)
ISM.D4 = na.omit(ISM.D4)
ADM.D0 = na.omit(ADM.D0)
ADM.D4 = na.omit(ADM.D4)
length(iPS$pbid)

isoforms = c(length(iPS$pbid), length(S1$pbid), length(S2$pbid), length(S3$pbid), length(ISM.D0$pbid), length(ISM.D4$pbid), length(ADM.D0$pbid), length(ADM.D4$pbid))
Stage = factor(Stage, levels = c("iPS","S1","S2","S3", "ISM.D0", "ISM.D4", "ADM.D0","ADM.D4"))
readsDistribut = data.frame(Stage, isoforms)
ggplot(readsDistribut, aes(x = Stage, y = isoforms, fill = Stage, label = reads)) + geom_bar(stat = "identity") +
  geom_text(size = 4,position = position_stack(vjust = 0.5))

gene1 = data.frame(Stage, isoforms_1)
View(gene1)
isoforms_1 = c(3996, 3604, 3758, 3703, 3497, 3757, 3788, 3679)
isoforms_over1=c(4869, 5516, 3364, 3666, 3179, 3277, 5350, 4846)
ggplot(gene1, aes(x = Stage, y = isoforms_1, fill = Stage, label = isoforms_1)) + geom_bar(stat = "identity") +
  geom_text(size = 4,position = position_stack(vjust = 0.5)) + labs(y = "gene numbers", title = "genes with 1 isoform")

gene2 = data.frame(Stage, isoforms_over1)
ggplot(gene2, aes(x = Stage, y = isoforms_over1, fill = Stage, label = isoforms_over1)) + geom_bar(stat = "identity") +
  geom_text(size = 4,position = position_stack(vjust = 0.5)) + labs(y = "gene numbers", title = "genes with more than 1 isoform")

stage = na.omit(iPS)
length(iPS$pbid)
freq = table(stage$refgene)
freq = as.data.frame(freq)
table(freq$Freq)
over1 = subset(freq, Freq>5)
sum(over1$Freq)
ggplot(freq, aes(x = Freq)) + geom_histogram(binwidth = 0.5, position = "identity", fill = "skyblue") + labs(x = "isoforms numbers", title = "ADM.D4") + xlim(0,50)
View(freq)

####intersect####

countADM.D4 = table(ADM.D4$refgene)
countISM.D4 = table(ISM.D4$refgene)
countS3 = table(S3$refgene)
length(countS3)
over1 = function(stage){
  i = 0
  over1gene = list()
  while(i<length(stage)){
    i = i+1
    if (stage[i]>1){
      over1gene = c(over1gene,names(stage[i]))
    }
  }
  return(over1gene)
}

ADM.D4over = over1(countADM.D4)
ADM.D4over = unlist(ADM.D4over)
ISM.D4over = over1(countISM.D4)
ISM.D4over = unlist(ISM.D4over)
S3over = over1(countS3)
S3over = unlist(S3over)
gene.ADM.D0 = unique(ADM.D0$refgene)
gene.ADM.D4 = unique(ADM.D4$refgene)
gene.ISM.D0 = unique(ISM.D0$refgene)
gene.ISM.D4 = unique(ISM.D4$refgene)
gene.S1 = unique(S1$refgene)
gene.S2 = unique(S2$refgene)
gene.S3 = unique(S3$refgene)
gene.iPS = unique(iPS$refgene)
gene.S2 = intersect(gene.S2,gene.S2)
genenumbers = c(max(iPS$pbgene), max(ADM.D0$pbgene), max(ADM.D4$pbgene), max(ISM.D0$pbgene), max(ISM.D4$pbgene), max(S1$pbgene), max(S2$pbgene), max(S3$pbgene))
gene.ISM.D0less3 = unique(ISM.D0less3$refgene)
gene.ADM.D0less3 = unique(ADM.D0less3$refgene)
Stage = c("ips", "ADM.D0", "ADM.D4", "ISM.D0","ISM.D4","S1","S2","S3", "ips", "ADM.D0", "ADM.D4", "ISM.D0","ISM.D4","S1","S2","S3","ips", "ADM.D0", "ADM.D4", "ISM.D0","ISM.D4","S1","S2","S3","ips", "ADM.D0", "ADM.D4", "ISM.D0","ISM.D4","S1","S2","S3","ips", "ADM.D0", "ADM.D4", "ISM.D0","ISM.D4","S1","S2","S3")
Stage = c("iPS","S1","S2","S3", "ISM.D0", "ISM.D4", "ADM.D0","ADM.D4")
numbers_stage = data.frame(stage, genenumbers)
ggplot(numbers_stage, aes(fill = stage, y = genenumbers, x = stage, label = genenumbers)) + geom_bar(stat = "identity") +
  geom_text(size = 4,position = position_stack(vjust = 0.5))

### find similarity
similarity = function(stage1,stage2){
  stage1 = unique(stage1)[(unique(stage1))]
  stage2 = unique(stage2)[(unique(stage2))]
  i = 0
  for (gene in stage1){
    if (gene %in% stage2){
      i = i+1
    }
  }
  return(i)
}

a = intersect(gene.iPS,gene.ADM.D0)
b = intersect(gene.S1, gene.ADM.D4)
c = intersect(gene.S2,gene.S3)
d = intersect(gene.ISM.D0,gene.ISM.D4)
e = intersect(a,b)
f = intersect(c,d)
g = intersect(e,f)
v = as.matrix(gene.iPS,gene.S1, gene.ADM.D0,gene.ADM.D4)
write.csv(g,"householdgene.txt", sep = "\r")
stage = c(iPS$refgene, ADM.D0$refgene, ADM.D4$refgene, ISM.D0$refgene, ISM.D4$refgene, S1$refgene, S2$refgene, S3$refgene)
m = 1
while (m<8){
  n = m+1
  if (n <= 8){
    sim = similarity(stage[m], stage[n])
    print(sim)
    n = n+1
  }
  m = m+1
}
stage[1]
similarity(S2,S3)
similarity("iPS")
### venn diagram
grid.newpage()
venn.plot <- draw.quad.venn(
  area1 = as.numeric(length(gene.iPS)),
  area2 = as.numeric(length(gene.S1)),
  area3 = as.numeric(length(gene.ADM.D0)),
  area4 = as.numeric(length(gene.ADM.D4)),
  n12 = as.numeric(length(intersect(gene.iPS,gene.S1))),
  n13 = as.numeric(length(intersect(gene.iPS, gene.ADM.D0))),
  n14 = as.numeric(length(intersect(gene.iPS,gene.ADM.D4))),
  n23 = as.numeric(length(intersect(gene.S1, gene.ADM.D0))),
  n24 = as.numeric(length(intersect(gene.S1,gene.ADM.D4))),
  n34 = as.numeric(length(intersect(gene.ADM.D0,gene.ADM.D4))),
  n123 = as.numeric(length(intersect(c,gene.ADM.D0))),
  n124 = as.numeric(length(intersect(c,gene.ADM.D4))),
  n134 = as.numeric(length(intersect(a,gene.ADM.D4))),
  n234 = as.numeric(length(intersect(d,gene.ADM.D4))),
  n1234 = as.numeric(length(intersect(c,f))),
  category = c("iPS", "S1", "ADM.D0", "ADM.D4"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
);


venn.diagram(
  x = list(gene.S2 , gene.ISM.D0less3 , gene.ADM.D0less3),
  category.names = c("S2" , "ISM.D0" , "ADM.D0"),
  filename = 'S2_ADM_ISM_less.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('skyblue', 'pink1', 'mediumorchid'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
##
numbers = c(2177,1774,1640,1950,1640,2513,3091,1966)
alt5 = c(161,224,167,228,267,384,3598,267)
alt3 = c(193,261,194,258,194,434,3756,276)
ES = c(33,37,63,58,63,65,285,36)
IR = c(13,11,21,10,21,28,68,14)
total = c(400,533,645,554,645,911,7707,593)
kinds = c(rep(c("alt5", "alt3", "ES", "IR", "total"), each = 8))
numbers = c(alt5, alt3, ES, IR, total)
stats = data.frame(Stage,numbers, kinds)
stats = data.frame(stage, numbers)
ggplot(data = stats,aes(x = stage, y = numbers)) +geom_bar(stat = "identity")
ggplot(data = stats,aes(x = Stage, y = numbers,fill = kinds)) + geom_bar(stat = "identity", position = position_dodge())
View(stats)

novel_gene_numbers =c(25,19,14,23,13,18,32,12)
confident = c()
genes = data.frame(stage,novel_gene_numbers)
ggplot(data = genes, aes(x = stage, y = novel_gene_numbers)) + geom_bar(stat = "identity")



### common genes