library(dada2)
library(data.table)

# fastq files reading
path <- "./inputDADA2"

fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names=T))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_cut"), `[`, 1)

# Quality visualisation
plotQualityProfile(fnFs[100:101])
plotQualityProfile(fnRs[114:115])

# Filtering and cutting reads according to their quality
filtFs <- file.path(".", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(".", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# truncLen values for each marker
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(170,150), maxN=0, truncQ=0, multithread=T, rm.phix=F)

# Error modelling
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

# Visualize error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Data dereplication
derepF <- derepFastq(filtFs)
derepR <- derepFastq(filtRs)

# Denoising
dadaF <- dada(derepF, err=errF, multithread=T, pool=F)
dadaR <- dada(derepR, err=errR, multithread=T, pool=F)

# Test for minOverlap value
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = 1, returnRejects = T)
overlap=c()
for (sam in sample.names) {
  overlap=c(overlap, mergers[[sam]]$nmatch)
}
plot(sort(overlap))

# Merge pairs
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, maxMismatch=0, minOverlap = 50)

# Length selection
seqtab <- makeSequenceTable(mergers)
table(nchar(colnames(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 200:210]

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=T, verbose=T)
sum(seqtab.nochim)/sum(seqtab2)

# Exporting data files
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "size", "nonchim")
rownames(track) <- sample.names
write.table(file="tab_track_dada2.csv", track, col.names=T, row.names=T, sep="\t", dec=".", quote=F)

ASVs=data.frame(Code=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_"), Sequence=getSequences(seqtab.nochim))
colnames(seqtab.nochim)=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_")
write.table(file="List_ASVs.fasta", ASVs, col.names=F, row.names=F, sep="\n", dec=".", quote=F)
write.table(file="tab_distri_ASVs.csv", seqtab.nochim, col.names=T, row.names=T, sep="\t", dec=".", quote=F)
