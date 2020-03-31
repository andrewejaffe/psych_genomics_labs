#########################
## explore GWAS output ##
#########################

### first run `plink_gwas.sh`
### through the frequency step

### install once
# install.packages("tidyverse")
# install.packages("qqman")

## initiate libraries
library(qqman)
library(tidyverse)

## read in genotype data
genotypes  = read_delim("example/wgas2.traw",delim = "\t")

## split up
snp = as.data.frame(genotypes[,-(1:6)])
map = as.data.frame(genotypes[,1:6])
rownames(snp) = map$SNP
rownames(map) = map$SNP

#################
## read in fam ##
# https://www.cog-genomics.org/plink2/formats#fam
# Sex code ('1' = male, '2' = female, '0' = unknown)
pheno = read_delim("example/wgas2.fam",delim=" ",
          col_names=c("FID", "IID", "FatherID","MotherID","Sex", "Y"))

## and MDS...this is formatted weird
mds = read_delim("example/wgas2.mds", delim=" ")
mds = read.table("example/wgas2.mds",as.is=TRUE,header=TRUE)

## join the two 
pheno = left_join(pheno, mds, by = c("FID", "IID"))

## make outcome a factor
pheno$Dx = factor(ifelse(pheno$Y == "1", "Control", "Case"),
                  levels = c("Control", "Case"))

## plots
qplot(x=C1,y=C2,data=pheno,col=Dx)
qplot(x=C1,y=C2,data=pheno,facets=~Dx)
qplot(x=C3,y=C4,data=pheno,col=Dx)

#####################
## read in results ##
#####################

lr= read.table("example/assoc1.assoc.logistic",as.is=TRUE,header=TRUE)

## diagnosis effect
outDx = lr %>% filter(TEST == "ADD")
rownames(outDx) = outDx$SNP
identical(rownames(outDx), rownames(map)) # TRUE

## how many missing
outDx %>% summarise(numNA = sum(is.na(OR)))
table(is.na(outDx$OR))

## by chr
outDx %>% group_by(CHR) %>% summarise(numNA = sum(is.na(OR)))

## check a few 3x2 tables
ii = which(is.na(outDx$P))[1]
table(snp[ii,], pheno$Dx)
## MAF effect?
maf = read.table("example/wgas2.frq", as.is=TRUE, header=TRUE)
map = left_join(map, maf)
## A1/A2 and COUNTED/ALT are two different ways
## plink codes alleles
identical(map$A1, map$COUNTED) # TRUE

## merge into results
outDx = map %>% select(SNP,MAF) %>% left_join(outDx, .)

## check MAFs
outDx %>% group_by(is.na(OR)) %>% summarise(mean(MAF))
outDx %>% group_by(is.na(OR)) %>% summarise(max(MAF))

####################
## manhattan plot ##
####################
plot(-log10(outDx$P))

outDx %>% filter(!is.na(P)) %>% manhattan()
outDx %>% filter(!is.na(P)) %>% manhattan(ylim=c(0,9))

## p-value vs maf
qplot(x = MAF, y=-log10(P), data=outDx)
qplot(x = cut(MAF,breaks = seq(0,0.5,by=0.05)), 
        y=-log10(P), data=outDx, geom="boxplot",
      xlab= "MAF category")

####################
#### QQ PLot ######

qq(outDx$P)

##########################
## component 1 effect ####
outC1 = lr %>% filter(TEST == "C1")
outC1 %>% filter(!is.na(P)) %>% manhattan(ylim=c(0,10))
## some kind of simulated dataset

qq(outC1$P)
