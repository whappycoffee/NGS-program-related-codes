novelgenes = "tapis_out/novelGenes.csv"
outdir = ["iPS/", "S1/", "S2/", "S3/", "ISM.D0/", "ISM.D4/", "ADM.D0/", "ADM.D4/"]

def writepos(stage):
	f = open(stage+novelgenes, "rb")
	lines = f.readlines()
	location = open(stage + "genepos.txt", "w+")
	for line in lines:
		pos = line.split()
		location.write(pos[0] + ":" + pos[1] + "-" + pos[2] + "\n")
	

for stage in outdir:
	writepos(stage)
