import subprocess
stages = ["iPS", "S1", "S2", "S3", "ISM.D0", "ISM.D4", "ADM.D0", "ADM.D4"]
cmd = "samtools view -c -F 260 "
def count(stage):
	if stage == "ADM.D4":
		cmd1 = cmd + stage + "/ToFU/hq_isoforms.fasta.sam"
		cmd2 = cmd + stage + "/ToFU/iso.collapsed.rep.fa.sam"
	else:
		cmd1 = cmd + stage + "/ToFU/hq_isoforms.fastq.sam"
		cmd2 = cmd + stage + "/ToFU/iso.collapsed.rep.fq.sam"
	print stage + "\r"
	status = subprocess.call(cmd1,shell=1)
	status = subprocess.call(cmd2,shell=1)

for stage in stages:
	count(stage)