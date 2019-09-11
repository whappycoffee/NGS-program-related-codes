library(VennDiagram)
ADM.D0 = read.csv("/Users/wuyibo/Documents/Research/results/ADM.D0/AS_summary.csv", header = TRUE, sep = ",")
ADM.D4 = read.csv("/Users/wuyibo/Documents/Research/results/ADM.D4/AS_summary.csv", header = TRUE, sep = ",")
iPS = read.csv("/Users/wuyibo/Documents/Research/results/iPS/AS_summary.csv", header = TRUE, sep = ",")
ISM.D0 = read.csv("/Users/wuyibo/Documents/Research/results/ISM.D0/AS_summary.csv", header = TRUE, sep = ",")
ISM.D4 = read.csv("/Users/wuyibo/Documents/Research/results/ISM.D4/AS_summary.csv", header = TRUE, sep = ",")
S1 = read.csv("/Users/wuyibo/Documents/Research/results/S1/AS_summary.csv", header = TRUE, sep = ",")
S2 = read.csv("/Users/wuyibo/Documents/Research/results/S2/AS_summary.csv", header = TRUE, sep = ",")
S3 = read.csv("/Users/wuyibo/Documents/Research/results/S3/AS_summary.csv", header = TRUE, sep = ",")

## all genes

ISM.D4 = read.csv("../ISM.D4/tapis_out/polyA_summary.csv", header = TRUE, sep = "\t")
ADM.D4 = read.csv("../ADM.D4/tapis_out/polyA_summary.csv", header = TRUE, sep = "\t")
ADM = rbind(ADM.D0,ADM.D4)
ISM = rbind(ISM.D0,ISM.D4)
length(S2$Gene)
S25 = subset(S2, S2$Total >= 5)
ADM.D45 = subset(ADM.D4, ADM.D4$Total >= 10)

ISM.D45 = subset(ISM.D4, ISM.D4$Total >=10)
length(ISM5$Gene)
length(unique(ISM5$Gene))
length(unique(ADM5$Gene))
venn.diagram(
  x = list(as.character(iPS$Gene) , as.character(ADM.D4$gene) , as.character(ISM.D4$gene)),
  category.names = c("iPS" , "ADM.D4" , "ISM.D4"),
  filename = 'iPS_ADMD4_ISMD4.png',
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

draw.pairwise.venn(1282, 1465, 549, category = c("ISM.D4", "ADM.D4"), lty = rep("blank", 
                                                                                   2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                        0), cat.dist = rep(0.025, 2), scaled = FALSE)
