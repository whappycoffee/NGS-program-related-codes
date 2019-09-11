/Users/wuyibo/src/gffread-0.10.8.OSX_x86_64/./gffread iso.collapsed.gff -T -o iso.collapsed.gtf
gene_model_to_splicegraph.py -m iso.collapsed.gtf -a -S -d predicted
find predicted -name "*.gff" > predicted.lis
genewise_statistics.py predicted.lis -C > AS_summary.csv