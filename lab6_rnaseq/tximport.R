## tximport

### package info
## https://bioconductor.org/packages/release/bioc/html/tximport.html
## https://f1000research.com/articles/4-1521/v2

### vignette
## https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
## https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.R

## ------------------------------------------------------------------------
library(tximportData)
dir <- system.file("extdata", package="tximportData")
list.files(dir)

## ------------------------------------------------------------------------
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
all(file.exists(files))

## ------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
# https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # make variable shorter
columns(txdb)
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

## ------------------------------------------------------------------------
library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)
## you can also do this with biomaRt
## https://bioconductor.org/packages/release/bioc/html/biomaRt.html
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene2txsym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","ensembl_transcript_id"), 
                   mart=ensembl)


## ------------------------------------------------------------------------
library(tximport)
### this may require more memory that the 1Gb on RStudio Cloud
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi$counts)

## ------------------------------------------------------------------------
## samlon, transcript-level
# txi.tx <- tximport(files, type="salmon", txOut=TRUE)
# txi.sum <- summarizeToGene(txi.tx, tx2gene)
# all.equal(txi$counts, txi.sum$counts)

## ------------------------------------------------------------------------
## these are for sailfish files
# tx2knownGene <- read_csv(file.path(dir, "tx2gene.csv"))
# files <- file.path(dir,"sailfish", samples$run, "quant.sf")
# names(files) <- paste0("sample",1:6)
# txi.sailfish <- tximport(files, type="sailfish", tx2gene=tx2knownGene)
# head(txi.sailfish$counts)


## ------------------------------------------------------------------------
## RSEM, gene level
# files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results.gz"))
# names(files) <- paste0("sample",1:6)
# txi.rsem <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)
# head(txi.rsem$counts)

## ---- results="hide", messages=FALSE-------------------------------------
library(edgeR)
library(csaw)

## ------------------------------------------------------------------------
cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length,
# adjusted to avoid changing the magnitude of the counts.
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- cts / normMat

# Computing effective library sizes from scaled counts,
# to account for composition biases between samples.
library(edgeR)
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors,
# and calculating offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
y <- y[keep,]
# y is now ready for estimate dispersion functions
# see edgeR User's Guide

## ------------------------------------------------------------------------
se <- SummarizedExperiment(assays=list(counts=y$counts, offset=y$offset))
se$totals <- y$samples$lib.size
library(csaw)
cpms <- calculateCPM(se, use.offsets=TRUE, log=FALSE)

## ---- results="hide", messages=FALSE-------------------------------------
library(DESeq2)

## ------------------------------------------------------------------------
sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=3)))
rownames(sampleTable) <- colnames(txi$counts)

## ------------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
# dds is now ready for DESeq()
# see DESeq2 vignette

## ------------------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
txi <- tximport(files, type="salmon",
                tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)
# filtering
keep <- filterByExpr(y)
y <- y[keep,]
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data=sampleTable)
v <- voom(y, design)
# v is now ready for lmFit()
# see limma User's Guide

## ------------------------------------------------------------------------
sessionInfo()
