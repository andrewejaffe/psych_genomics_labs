## using already processed data with recount

# https://jhubiostatistics.shinyapps.io/recount/
# https://f1000research.com/articles/6-1558
# https://www.nature.com/articles/nbt.3838

## learn more about [Ranged]SummarizedExperiments
# https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
library(SummarizedExperiment)

# data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32465
f = "http://duffel.rail.bio/recount/v2/SRP009615/rse_gene.Rdata"
utils::download.file(f, destfile="SRP009615_rse_gene.RData", mode="wb") 

## load it
load("SRP009615_rse_gene.RData",verbose = TRUE)
rse_gene

## add phenotype data
library(readr)
sample_info = read_csv("SRP009615_sample_info.csv")

## add sample info
colData(rse_gene)$group <- sample_info$group
colData(rse_gene)$gene_target <- sample_info$gene_target

### scale counts
## from the recount package
source("scale_counts.R")
rse <- scale_counts(rse_gene)

## do differential expression
library('DESeq2')

## standard DESeq2 analysis
dds <- DESeqDataSet(rse, ~ gene_target + group)
dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target, fitType = 'local')
res <- results(dds)

## make plots
plotMA(res, main="DESeq2 results for SRP009615")

################
# with limma
library(edgeR)

mod = model.matrix(~ gene_target + group, data=colData(rse_gene))

## filter expression
keepIndex = rowSums(assays(rse)$counts) > 30
rse_filter = rse[keepIndex,]

## fix symbols
rowData(rse_filter)$symbol = sapply(rowData(rse_filter)$symbol,"[",1)

##### DGE ######
dge = DGEList(counts = assays(rse_filter)$counts, 
              genes = rowData(rse_filter))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)

## do analysis
fitGene = lmFit(vGene)

## top table
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
                   p.value = 1,number=nrow(rse), sort="none")
sum(outGene$adj.P.Val < 0.05)

outGene[outGene$adj.P.Val < 0.05,]

