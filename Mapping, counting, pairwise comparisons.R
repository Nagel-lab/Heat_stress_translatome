
library(systemPipeR)


## Read mapping with `HISAT2`
#----------------------------

args <- systemArgs(sysma="hisat2.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
moduleload(modules(args))
system("hisat2-build ./general/tair10.fasta ./general/tair10.fasta")
resources <- list(walltime=120, ntasks=1, ncpus=5, memory=10240) 
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs=53, template = "batchtools.slurm.tmpl", runid="01", resourceList=resources)


## Read and alignment stats
#---------------------------

read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "Results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


# Read quantification
#---------------------

## Read counting with `summarizeOverlaps` in parallel mode using multiple cores

library("GenomicFeatures"); library(BiocParallel)
txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
saveDb(txdb, file="./data/tair10.sqlite")
txdb <- loadDb("./tair10.sqlite")
(align <- readGAlignments(outpaths(args)[1])) # Demonstrates how to read bam file into R
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=2); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=TRUE)) 
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "Results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "Results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


## Sample-wise correlation analysis
#----------------------------------

library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
countDF <- as.matrix(read.table("./Results/countDFeByg.xls"))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("Results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()


## Export rlog values
#--------------------
rlog <- assay(rlog(dds))
write.table(rlog, file = "./Results/rlog_values.txt", sep = "\t")


## Analysis of DEGs
#-------------------
countDF <- read.delim("Results/countDFeByg.xls", row.names=1, check.names=FALSE) 
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)


## Add gene descriptions
#----------------------

library("biomaRt")
m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
desc <- getBM(attributes=c("tair_locus", "description"), mart=m)
desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
degseqDF <- data.frame(degseqDF, Desc=descv[rownames(degseqDF)], check.names=FALSE)
write.table(degseqDF, "./Results/DESeq2_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)


## Plot DEG results
#------------------
degseqDF <- read.delim("Results/DESeq2_allcomp.xls", row.names=1, check.names=FALSE) 
pdf("Results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list$Summary, "./Results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)

