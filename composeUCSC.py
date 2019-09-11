import subprocess
Dir = "/Users/wuyibo/Documents/Research/results/"
stages = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]
def command(stage):
	#cmd1 = "cd " + Dir + stage
	cmd2 = "samtools sort -o " + stage + "/" + stage + ".sorted.bam " + stage + "/" + stage + ".bam"
	cmd3 = "samtools index " + stage + "/" + stage + ".sorted.bam"
	cmd4 = "samtools view " + stage + "/" + stage + ".sorted.bam X -b > " + stage + "/" + stage + ".X.bam"
	cmd5 = "samtools index " + stage + "/" + stage + ".X.bam"
	#subprocess.call(cmd1,shell=1)
	subprocess.call(cmd2,shell=1)
	subprocess.call(cmd3,shell=1)
	subprocess.call(cmd4,shell=1)
	subprocess.call(cmd5,shell=1)
for stage in stages:
	command(stage)