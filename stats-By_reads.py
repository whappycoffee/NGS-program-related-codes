import pandas as pd 
import numpy as np 
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from matplotlib import style
import operator
style.use("ggplot")
#stage = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]
Dir = "/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt"

def get_counts(stage):
	df = pd.read_csv(stage+Dir, sep = "\t")
	df.dropna(inplace = True)
	counts = Counter(df["refgene"])
	return counts
S1 = get_counts("S1")
S2 = get_counts("S2")
S3 = get_counts("S3")
ADMD0 = get_counts("ADM.D0")
ADMD4 = get_counts("ADM.D4")
ISMD0 = get_counts("ISM.D0")
ISMD4 = get_counts("ISM.D4")
iPS = get_counts("iPS")

#colors = ['#E69F00', '#56B4E9']
#names = ["iPS", "S1"]
#def genenumbers(stage):
#	iso1 = 0
#	iso2 = 0
#	isonumbers = 0
#	for gene in stage:
#		isonumbers+=stage[gene]
#		if stage[gene] == 1:
#			iso1+=1
#		elif stage[gene] > 1:
#			iso2+=1
#	return iso1, iso2, isonumbers#

#stages = [iPS, S1, S2, S3, ISMD0, ISMD4, ADMD0, ADMD4]
#for stage in stages:
#	iso1, iso2,isonumbers = genenumbers(stage)
#	print "1 \t", iso1, "\t 2 \t", iso2, "\t total \t", isonumbers
def common_genes(stage1,stage2):
	genenumbers = 0
	isoforms = []
	isoforms2 = []
	total = 0
	for gene in stage1:
		if gene in stage2:
			total+=1
		if gene in stage2 and abs(stage1[gene] - stage2[gene])<=5:
			genenumbers+=1
			isoforms.append(stage1[gene])
			isoforms2.append(stage2[gene])
	SUM = sum(Counter(isoforms).values())
	less5 = Counter(isoforms)[1] + Counter(isoforms)[2]+Counter(isoforms)[3]+Counter(isoforms)[4]+Counter(isoforms)[5]
	return genenumbers
	#return (genenumbers*100.00)/(len(stage1) + len(stage2) - total)

#def merge(stage1, stage2):
#	newdata = stage1
#	for gene in stage2:
#		if gene in stage1:
#			newdata.update({gene:max([stage1[gene], stage2[gene]])})
#		if gene not in stage1:
#			newdata.update({gene:stage2[gene]})
#	return newdata
#ISM = merge(ISMD0,ISMD4)
#ADM = merge(ADMD0,ADMD4)#

def sortdic(stage):
	sort_x = sorted(stage.items(),key=operator.itemgetter(1))
	return sort_x
def commonTopgene(stage1,stage2,stage3):
	G = []
	for gene in stage1:
		if gene in stage2 and gene in stage3:
			G.append(gene)
	return G

print commonTopgene(sortdic(S3)[-10:], sortdic(ISMD4)[-10:], sortdic(ADMD4)[-10:])
print commonTopgene(sortdic(S2)[-10:], sortdic(ISMD0)[-10:], sortdic(ADMD0)[-10:])


#print sortdic(iPS)[-10:]
#print sortdic(S1)[-10:]
#print sortdic(S2)[-10:]
print sortdic(S3)[-10:]
#print sortdic(ADMD0)[-10:]
print sortdic(ADMD4)[-10:]
#print sortdic(ISMD0)[-10:]
print sortdic(ISMD4)[-10:]

#print "ISMD0 AMDD0", common_genes(ISMD0,ADMD0)
#print "ISMD0 AMDD4", common_genes(ISMD0,ADMD4)
#print "ISMD4 AMDD0", common_genes(ISMD4,ADMD0)
#print "ISMD4 AMDD4", common_genes(ISMD4,ADMD4)
#print "ISM S2", common_genes(ISMD0,S2)
#print "ISM S2", common_genes(ISMD4,S2)

##a = [42.48,42.11,42.35,42.14,51.59,44.32,43.67,51.70,41.86,42.51,48.21,42.05,44.64,42.20,44.81,47.21,42.30,46.25,49.35,45.56,49.65,47.79,44.31,48.32,50.68,47.33,47.07,47.51]
#print sorted(a)
#iso1, iso2 = common_genes(iPS, S3)
#def less5proportion(iso1, iso2):
#	SUM = sum(Counter(iso1).values())
#	less5 = Counter(iso1)[1] + Counter(iso1)[2]+Counter(iso1)[3]+Counter(iso1)[4]+Counter(iso1)[5]
#	return (less5*100.00)/SUM

#print less5proportion(iso1,iso2)
#plt.hist(iso1_iPS,color = colors[1], edgecolor = 'black', bins = 22)
#plt.legend()
#plt.xlabel("isoform numbers")
#plt.ylabel("counts")
#plt.title("isoform number distribution of iPS-S1 common genes")
#plt.show()

#for gene in iPS:
#	if gene in ADMD0 and gene in ADMD4 and gene in ISMD0 and gene 
#print "S1, ADM.D4", ((common_genes(S1, ADMD4)*100.00)/(len(S1)+len(ADMD4)))


#print "S1, ADM.D4", common_genes(S1, ADMD4)
#print "S1, ADM.D0", common_genes(S1, ADMD0)
#print "S1, ISM.D0", common_genes(S1, ISMD0)
#print "S1, ISM.D4", common_genes(S1, ISMD4)
#print "S1, iPS", common_genes(S1, iPS)
#print "S1, S2", common_genes(S1, S2)
#print "S1, S3", common_genes(S1, S3)
#print "S2, ADM.D4", common_genes(S2, ADMD4)
#print "S2, ADM.D0", common_genes(S2, ADMD0)
#print "S2, ISM.D0", common_genes(S2, ISMD0)
#print "S2, ISM.D4", common_genes(S2, ISMD4)
#print "S2, iPS", common_genes(S2, iPS)
#print "S2, S3", common_genes(S2, S3)
#print "S3, ADM.D4", common_genes(S3, ADMD4)
#print "S3, ADM.D0", common_genes(S3, ADMD0)
#print "S3, ISM.D0", common_genes(S3, ISMD0)
#print "S3, ISM.D4", common_genes(S3, ISMD4)
#print "S3, iPS", common_genes(S3, iPS)
#print "ADM.D0, ADM.D4", common_genes(ADMD0, ADMD4)
#print "ADM.D0, ISM.D0", common_genes(ADMD0, ISMD0)
#print "ADM.D0, ISM.D4", common_genes(ADMD0, ISMD4)
#print "ADM.D0, iPS", common_genes(ADMD0, iPS)
#print "ADM.D4, ISM.D0", common_genes(ADMD4, ISMD0)
#print "ADM.D4, ISM.D4", common_genes(ADMD4, ISMD4)
#print "ADM.D4, iPS", common_genes(ADMD4, iPS)
#print "ISM.D0, ISM.D4", common_genes(ISMD0, ISMD4)
#print "ISM.D0, iPS", common_genes(ISMD0, iPS)
#print "ISM.D4, iPS", common_genes(ISMD4, iPS)


#df = pd.read_csv("iPS/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt", sep = "\t")
#df.dropna(inplace = True)
#counts = Counter(df["refgene"])
#print counts
# print Counter(df["refgene"]).items()[1:5]
# table = np.array([[0,0]])
# for key in counts:
# 	table = np.append(table,[[key, counts[key]]],axis = 0)
# table = np.delete(table, 0, axis = 0)
# newdf = pd.DataFrame(table, columns = ["genes", "counts"])
# print newdf.head()
#sns.distplot(freq, kde = False, rug = True)
#print freq
#plt.show()
#genes = pd.Series([tuple(i) for i in df["refgene"]])
#counts = genes.value_counts()
#print counts[0:5]





#print Counter(df["refgene"]).items()[1:5]