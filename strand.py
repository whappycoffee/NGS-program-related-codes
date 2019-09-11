import subprocess
Dir = "/Users/wuyibo/Documents/Research/results/"
stages = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]
def command(stage):
	#cmd1 = "cd " + Dir + stage
	cmd2 = "samtools view -c -F 16 " + stage + "/" + stage + ".sorted.bam X:31,115,794-33,357,558"
	cmd3 = "samtools view -c -f 16 " + stage + "/" + stage + ".sorted.bam X:31,115,794-33,357,558"
	print stage
	#subprocess.call(cmd1,shell=1)
	subprocess.call(cmd2,shell=1)
	subprocess.call(cmd3,shell=1)
	#subprocess.call(cmd4,shell=1)
	#subprocess.call(cmd5,shell=1)
for stage in stages:
	command(stage)
