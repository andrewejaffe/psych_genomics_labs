##############################
# Intro to Bioconductor lab ##
##############################

## installing bioconductor packages uses
##  different syntax than base R packages
##  read more: http://bioconductor.org/install/
## this has also changed relatively recently
## theres a ton of courses here: http://bioconductor.org/help/course-materials/

## you need the BiocManager CRAN package now
if (!requireNamespace("BiocManager"))  install.packages("BiocManager")

## more here in their vignette:
## https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html
BiocManager::install(c("GenomicRanges","Biostrings"))
##this installs a bunch of packages and takes a bit of time

## we can check out a few of these course materials while that installs
# http://bioconductor.org/help/course-materials/2017/OSU/B1_Bioconductor_intro.html#bioconductor-sequencing-ecosystem 

### the rest is modeled after some of these courses
# http://bioconductor.org/help/course-materials/2017/OSU/B2_Common_Operations.html
# http://bioconductor.org/help/course-materials/2017/BioCAsia/A01-Bioc-Basics.html

############################
#### exploring sequences ###
############################

# https://bioconductor.org/packages/release/bioc/vignettes/Biostrings/inst/doc/BiostringsQuickOverview.pdf
library(Biostrings)

dna <- DNAStringSet( c("AAACTCTTG", "CCTTCAACA") )
dna

## various operations
reverse(dna)
complement(dna)
reverseComplement(dna)
translate(dna) # must be multiple of length 3

## extraction
subseq(dna, 1, 4)

## summaries
alphabetFrequency(dna)
dinucleotideFrequency(dna)

## counting
vcountPattern("AA", dna) # when multiple sequences
vmatchPattern("AA", dna) # when multiple sequences

##############################
## IntegerRanges / IRanges ###
##############################

# https://www.bioconductor.org/packages/release/bioc/html/IRanges.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf

### code chunk number 5: iranges-constructor
ir1 <- IRanges(start=1:10, width=10:1)
ir1
ir2 <- IRanges(start=1:10, end=11)
ir3 <- IRanges(end=11, width=10:1)
identical(ir1, ir2) && identical(ir1, ir3)
ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
              width=c(12, 6, 6, 15, 6, 2, 7))
ir

## common accessors 
start(ir)
end(ir)
width(ir)
ir[1:4]
ir[start(ir) <= 15]

## making a custom function
### code chunk number 11: plotRanges
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...) {
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

## plotting example
plotRanges(ir)
## compare to reduce
reduce(ir)
plotRanges(reduce(ir))
## compare to disjoin
disjoin(ir)
plotRanges(disjoin(ir))

## coverage
coverage(ir)
as.numeric(coverage(ir))

#############################
## GenomicRanges / GRanges ##
#############################

# https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html

library(GenomicRanges)

## ----example-GRanges-------------------------------------------------------
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr

## ----GRanges-location-accessors--------------------------------------------
seqnames(gr)
ranges(gr)
strand(gr)

## ----granges-accessor------------------------------------------------------
granges(gr)

## ----metadataAccess--------------------------------------------------------
mcols(gr)
mcols(gr)$score
gr$score

## ----setSeqLengths---------------------------------------------------------
seqlengths(gr) <- c(249250621, 243199373, 198022430)

## ----setSeqLengths2--------------------------------------------------------
seqlengths(gr)

## ----names-----------------------------------------------------------------
names(gr)
length(gr)
width(gr)

## ----splitAppendGRanges----------------------------------------------------
sp <- split(gr, rep(1:2, each=5))
sp
start(sp)

## ----other-----------------------------------------------------------------
rev(gr)
head(gr,n=2)
tail(gr,n=2)
window(gr, start=2,end=4)
gr[IRanges(start=c(2,7), end=c(3,9))]

## ----IRangesStuff----------------------------------------------------------
singles <- split(gr, names(gr))
g <- gr[1:3]
g <- append(g, singles[[10]])

##
start(g)
end(g)
width(g)
range(g)

## ----flank-----------------------------------------------------------------
flank(g, 10)

## ----flank2----------------------------------------------------------------
flank(g, 10, start=FALSE)

## ----shiftAndResize--------------------------------------------------------
shift(g, 5)
shift(g, -5)
resize(g, 30)

## ----reduce----------------------------------------------------------------
reduce(g)

## ----gaps------------------------------------------------------------------
gaps(g)

## ----disjoin---------------------------------------------------------------
disjoin(g)

## ----coverage--------------------------------------------------------------
coverage(g)

## ----intervals1------------------------------------------------------------
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
setdiff(g, g2)

## ----example-overlaps---------------------------------------------------
gr1 <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(103, 106),
  strand = "+",
  score = 5L, GC = 0.45)

## ----countOL---------------------------------------------------------------
countOverlaps(gr, gr1)

## ----subsetByOverlaps------------------------------------------------------
subsetByOverlaps(gr,gr1)

## ----select-first----------------------------------------------------------
findOverlaps(gr, gr1, select="first")
findOverlaps(gr1, gr, select="first")
findOverlaps(gr, gr1)


