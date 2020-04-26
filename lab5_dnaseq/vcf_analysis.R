if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

library(VariantAnnotation)

fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
vcf

## header
header(vcf)

## sample names
samples(header(vcf))

## geno details
geno(header(vcf))


###################################################
head(rowRanges(vcf), 3)


## details of sequence
ref(vcf)[1:5]
qual(vcf)[1:5]
alt(vcf)[1:5]

## genotypes themselves
geno(vcf)
sapply(geno(vcf), class)


geno(header(vcf))["DS",]

DS <-geno(vcf)$DS
dim(DS)
DS[1:3,]
