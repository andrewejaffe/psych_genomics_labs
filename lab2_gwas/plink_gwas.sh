########################
## example plink GWAS ##
########################

## Click terminal
cd /cloud/project/lab2_gwas


## based on https://atgu.mgh.harvard.edu/plinkseq/gwas.shtml
## and https://www.cog-genomics.org/plink/1.9/resources#teach


## get plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip
wget https://www.cog-genomics.org/static/bin/plink/teaching.zip
wget http://zzz.bwh.harvard.edu/plink/dist/example.zip

## unzip
unzip plink_linux_x86_64_20190304.zip 
unzip teaching.zip
unzip example.zip

## remove ssome stuff
rm teaching.zip
rm plink_linux_x86_64_20190304.zip
rm practical-1-slides.ppt

###########################
## useful linux commands ##
###########################

## change prompt
PS1='\W> '
## need to run every session

head toy.map 
ls # lists contents of current directory
wc -l toy.map # number of rows
ls -l 
cd ..
# executing files
plink --help # ?
./plink --help 

# PATH
echo $PATH # linux variables
export PATH="$PATH:/cloud/project/lab2_gwas"
plink --help # there we go

##########################
#### plink data types ####

# plain text (map/ped) 
head example/wgas1.map
wc -l example/wgas1.map
more example/wgas1.map
# space to advance, ctrl+c to exit

###################
## convery to binary (bim/fam/bed)
plink --file example/wgas1 --make-bed --out example/wgas2 

ls -l example

head example/wgas2.fam
head example/wgas2.bim
head example/wgas2.bed # doesnt work since binarys

########################
## get MDS components ##
## for quant ethnicity #

# first get independent SNPs
plink --bfile example/wgas2 --indep 100 10 1.25 --out example/wgas2 
# VIF = 1/(1-R2) so 1.25 = 1/(1-0.2) means R^2 = 0.2

# and then cluster and calculate components
# note this is different from the example text
# since its new in plink2
plink --bfile example/wgas2 --cluster --mds-plot 5 \
  --extract example/wgas2.prune.in --out example/wgas2 

### GWAS example with logistic regression
plink --bfile example/wgas2 --logistic \
       --covar example/wgas2.mds \
       --covar-name C1-C2 --out example/assoc1 

## their `awk` example
awk ' NR == 1 || ( $5 == "ADD" && $9 < 1e-4 ) ' example/assoc1.assoc.logistic 
# in the first row, check if column 5 is 'ADD' and
# if column 9 (the p-value) is less than 1e-4 

## their `grep` example
grep -w ADD example/assoc1.assoc.logistic | head -3
# match rows that have the string ADD

## we can also make a "transposed" file 
## that we can read into R later
plink --bfile example/wgas2 --recode A-transpose --out example/wgas2

## calculate MAF / minor allele frequency
plink --bfile example/wgas2 --freq --out example/wgas2 

################################
######### END HERE FOR MAIN ####
################################

## plinkseq examples ####

## add to path, otherwise you will need to type
## plinkseq-0.10/pseq every time you want to run
export PATH="$PATH:plinkseq-0.10"
pseq --help

## new project
pseq proj new-project 

## load data
pseq proj load-plink --file example/wgas2 --id gwas 

## add mds components
## this is just their code
awk ' BEGIN { printf "##mds1,Float,X,\"MDS component 1\"\n"; \
              printf "##mds2,Float,X,\"MDS component 2\"\n"; \
              printf "#ID\tmds1\tmds2\n" } \
       NR>1 { printf $1"\t"$4"\t"$5"\n" } ' example/wgas2.mds \
 > example/mds1.phe 
### now add 
pseq proj load-pheno --file example/mds1.phe 

## do logistic
pseq proj glm --phenotype phe1 --covar mds1 mds2 
