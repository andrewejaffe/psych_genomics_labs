## first install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi", version = "3.8")
BiocManager::install("limma", version = "3.8")
BiocManager::install("IlluminaHumanMethylation450kmanifest", version = "3.8")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", version = "3.8")
BiocManager::install("FlowSorted.DLPFC.450k", version = "3.8")
install.packages("tidyverse")

## libraries
library(minfi)
library(limma)
library(readr)
library(ggplot2)
# based on:
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html

## read in manifest
manifest= read_csv("pheno_n9_brainMeth.csv")
manifest = as.data.frame(manifest)
rownames(manifest) = manifest$BrNum

## check samples
table(manifest$Age < 0, manifest$Dx)
manifest$Group = ifelse(manifest$Age < 0, "Prenatal", 
	ifelse(manifest$Dx == "Schizo", "SCZD", "CONT"))
table(manifest$Group)

## make RGset
manifest$BasePath = paste0("idat/", manifest$Chip)
RGset = read.metharray(manifest$BasePath, verbose=TRUE)
colnames(RGset) = manifest$BrNum
metadata(RGset) = manifest

## check annotation
annotation(RGset)

## check probe info
probeMan = getManifest(RGset)
probeMan

head(getProbeInfo(probeMan))

## convert to methylset
Mset <- mapToGenome(RGset, mergeManifest  = TRUE) 
Mset

## get intensities for each channel
head(getMeth(Mset)[,1:3])
head(getUnmeth(Mset)[,1:3])

## convert to ratio set
Rset <- ratioConvert(Mset, what = "both", keepCN = TRUE)
Rset
beta <- getBeta(Rset)

##############
### QC #######
qc <- getQC(Mset)
plotQC(qc)
qcReport(RGset, sampNames = manifest$BrNum, 
	sampGroups = manifest$Group,pdf= "qcReport.pdf")

## and sex
predictedSex <- getSex(Mset, cutoff = -2)$predictedSex
head(predictedSex)
all.equal(manifest$Gender, predictedSex) # TRUE

###############
## normalize ##

## quantile
GRset.quantile <- preprocessQuantile(RGset, fixOutliers = TRUE,
  removeBadSamples = TRUE, badSampleCutoff = 10.5,
  quantileNormalize = TRUE, stratified = TRUE, 
  mergeManifest = TRUE, sex = NULL)
  
## noob
Mset.noob <- preprocessNoob(RGset)

## funnorm
GRset.funnorm <- preprocessFunnorm(RGset)

##################
## post process ##

## snps
snps <- getSnpInfo(GRset.quantile)
head(snps,10)

## add it
GRset.quantile <- addSnpInfo(GRset.quantile)

## you can also drop
GRset.quantile <- dropLociWithSnps(GRset.quantile, 
		snps=c("SBE","CpG"), maf=0)
GRset.quantile

## composition
# library(FlowSorted.DLPFC.450k)
# cellCounts <- estimateCellCounts(RGset, compositeCellType = "DLPFC")


###############
## analysis ###
###############

## add metadata
pData(GRset.quantile) = DataFrame(manifest)

## get meth
beta <- getBeta(GRset.quantile)
group  <- pData(GRset.quantile)$Group

#########
## first analysis, pre vs postnatal
devIndex = which(group %in% c("Prenatal", "CONT"))
betaDev = beta[,devIndex]
groupDev = factor(group[devIndex], levels = c("Prenatal", "CONT"),
	labels = c("Prenatal", "Postnatal"))

## with their wrapper
dmpDev <- dmpFinder(betaDev, pheno = groupDev,
	type = "categorical",shrinkVar=TRUE)
head(dmpDev)

## with limma
modDev = model.matrix(~groupDev)
fitDev = lmFit(betaDev, modDev)
ebDev = topTable(eBayes(fitDev), coef=2, n = nrow(betaDev), sort="none")
colnames(ebDev)[1:2] = c("diffMean", "AveMeth")

## look at results
qplot(x = ebDev$P.Value)
qplot(x = ebDev$P.Value, geom = "histogram",xlab="P-value")
qplot(x = ebDev$diffMean, geom = "histogram",xlab="diff in means")

#########
## second analysis, case vs adult control
dxIndex = which(group %in% c("SCZD", "CONT"))
betaDx = beta[,dxIndex]
groupDx = factor(group[dxIndex])

## with limma
modDx = model.matrix(~groupDx)
fitDx = lmFit(betaDx, modDx)
ebDx = topTable(eBayes(fitDx), coef=2, n = nrow(betaDx), sort="none")
colnames(ebDx)[1:2] = c("diffMean", "AveMeth")

## look at results
qplot(x = ebDx$P.Value, geom = "histogram",xlab="P-value")
qplot(x = ebDx$diffMean, geom = "histogram",xlab="diff in means")
