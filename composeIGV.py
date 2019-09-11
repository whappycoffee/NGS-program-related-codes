stages = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]

content = "sort strand" + "\r" + "expand Homo_sapiens.GRCh37.87.sorted.gtf" + "\r" + "snapshot"+ "\r"
def writebatch(stage):
	pos = open(stage + "/" + "genepos.txt", "rb")
	lines = pos.read().splitlines()[1:]
	IGV = open(stage + "/" + "IGV.txt", "w+")
	IGV.write("new" + "\r" + "genome hg19" + "\r" + "load /Users/wuyibo/Documents/Research/data/Homo_sapiens.GRCh37.87.sorted.gtf" + "\r" + "load /Users/wuyibo/Documents/Research/results/" + stage + "/" + stage + ".bam" + "\r" + "load /Users/wuyibo/Documents/Research/results/"+ stage +"/tapis_out/assembled.gtf" + "\r" + "snapshotDirectory /Users/wuyibo/Documents/Research/results/"+ stage +"/novelgenes/"+ "\r")
	for line in lines:
		numbers = line.split(":")
		loc = numbers[1].split("-")
		loci = numbers[0] + ":" + str(int(loc[0])-100) + "-" + str(int(loc[1])+100)
		IGV.write("goto " + loci + "\r" + content)

for stage in stages:
	writebatch(stage)
	
