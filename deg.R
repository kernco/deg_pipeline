require(DESeq2)
args = commandArgs(TRUE)
control = scan(snakemake@input[["control"]], sep=" ", what="")
treatment = scan(snakemake@input[["condition"]], sep=" ", what="")
sampleFiles = c(control, treatment)
sampleNames = sampleFiles
condition = c(rep("control", length(control)), rep("treatment", length(treatment)))
sampleTable = data.frame(sampleName=sampleNames[sampleNames != ""], fileName=sampleFiles[sampleFiles != ""], condition=condition)
dds = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=".",design=~condition)
dds = DESeq(dds)
res = results(dds)
resOrdered = res[order(res$padj),]
write.table(as.data.frame(resOrdered),quote=FALSE,sep='\t',file=snakemake@output)
