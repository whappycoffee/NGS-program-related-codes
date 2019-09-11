stages = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]
iPS = []
S1 = []
S2 = []
S3 = []
ADMD0 = []
ADMD4 = []
ISMD0 = []
ISMD4 = []
Stages = [iPS,S1,S2,S3,ADMD0,ADMD4,ISMD0,ISMD4]
f = open("ADM.D0/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
lines = f.readlines()
genenames = []
for line in lines:
	words = line.split("\t")
	genenames.append(words[3])

iPS = list(dict.fromkeys(genenames))

def getgenes(stage):
	if stage == "ADM.D4":
		f = open(stage+"/ToFU/iso.collapsed.rep.fa.sorted.sam.matchAnnot.txt.parsed.txt")
		lines = f.readlines()
		genenames = []
		for line in lines:
			words = line.split("\t")
			genenames.append(words[3])
	else:
		f = open(stage+"/ToFU/iso.collapsed.rep.fq.sorted.sam.matchAnnot.txt.parsed.txt")
		lines = f.readlines()
		genenames = []
		for line in lines:
			words = line.split("\t")
			genenames.append(words[3])
	return genenames
i = 0
for stage in stages:
	Stages[i] = getgenes(stage)
	i+=1
j = 0
for gene in iPS:
	if gene in ADMD0 and gene in ADMD4 and gene in ISMD0 and gene in ISMD4 and gene in S1 and gene in S2 and gene in S3:
		j+=1

print j



