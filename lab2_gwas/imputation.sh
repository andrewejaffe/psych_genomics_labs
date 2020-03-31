####################
# impute2 examples #
####################

cd /cloud/project/lab2_gwas

## get software
wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
tar xvfz impute_v2.3.2_x86_64_static.tgz
export PATH="$PATH:/cloud/project/lab2_gwas/impute_v2.3.2_x86_64_static"


###########################
## examples from impute2 ##
###########################

## change folders
cd impute_v2.3.2_x86_64_static

## imputation with 1 phased reference
## https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex1
./impute2 \
 -m ./Example/example.chr22.map \
 -h ./Example/example.chr22.1kG.haps \
 -l ./Example/example.chr22.1kG.legend \
 -g ./Example/example.chr22.study.gens \
 -strand_g ./Example/example.chr22.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 \
 -o ./Example/example.chr22.one.phased.impute2

# using shell variables, can save typing
E=/cloud/project/lab2_gwas/impute_v2.3.2_x86_64_static/Example/example.chr22
./impute2 -m $E.map -h $E.1kG.haps -l $E.1kG.legend \
 -g $E.study.gens -strand_g $E.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 -o $E.one.phased.shortcut.impute2

########
## imputation with 1 phased reference with pre-phasing
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex2

## first prephasing
./impute2 -prephase_g -m $E.map -g $E.study.gens \
 -int 20.4e6 20.5e6 -Ne 20000 -o $E.prephasing.impute2
## then impute
./impute2 -use_prephased_g -m $E.map -h $E.1kG.haps \
 -l $E.1kG.legend -known_haps_g $E.prephasing.impute2_haps \
 -strand_g $E.study.strand -int 20.4e6 20.5e6 \
 -Ne 20000 -o $E.one.prephased.impute2 -phase
 
#################
## imputation to chrX 
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex3

## this file doesn't exist the latest download
## but this is the gist of the code
## note the -chrX flag
X=/cloud/project/lab2_gwas/impute_v2.3.2_x86_64_static/Example/chrX/example.chrX
./impute2 -chrX \
 -m $X.map -h $X.reference.hap -l $X.reference.legend \
 -g $X.study.gen -sample_g $X.study.sample \
 -int 10.3e6 10.7e6 -Ne 20000 \
 -o $X.one.phased.impute2
 
###################
# one phased reference with variant filtering
## need to drop LOWCOV out of their filtering example
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex4
./impute2 \
 -filt_rules_l 'eur.maf<0.01' 'afr.maf<=0.05' \
 -m $E.map -h $E.1kG.haps -l $E.1kG.annot.legend \
 -g $E.study.gens -strand_g $E.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 -o $E.one.phased.filtered.impute2
 
######################
## unphased example
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex5

./impute2 \
 -m $E.map -g_ref $E.reference.gens \
 -strand_g_ref $E.reference.strand -g $E.study.gens \
 -strand_g $E.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 \
 -o $E.one.unphased.impute2
 
#####################
## imputation with two references
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex6
./impute2 \
 -m $E.map -h $E.1kG.haps $E.hm3.haps \
 -l $E.1kG.legend $E.hm3.legend -g $E.study.gens \
 -strand_g $E.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 \
 -o $E.two.phased.impute2
 
#######################
## just phasing #####
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex10
./impute2 -phase \
  -m $E.map -g $E.study.gens \
 -int 20.4e6 20.5e6 -Ne 20000 \
 -o $E.phasing.impute2

head -3 $E.phasing.impute2

## phasing with reference 
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex11
./impute2 -phase \
 -m $E.map -h $E.1kG.haps -l $E.1kG.legend -g $E.study.gens \
 -strand_g $E.study.strand \
 -int 20.4e6 20.5e6 -Ne 20000 \
 -o $E.phasing.ref.impute2



###########################
### strand alignment ######
# this is super important, see this:
# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#strand_options
## all of the above examples are already lined up

# shapeit is a good tool to check your data first
cd /cloud/project/lab_gwas
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz
tar xvfz shapeit.v2.r904.glibcv2.17.linux.tar.gz
export PATH="$PATH:/cloud/project/lab2_gwas/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin"

# back to impute2 folder
cd impute_v2.3.2_x86_64_static

## shapeit check
B=/cloud/project/lab2/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/example

## no reference
shapeit -check -G $B/gwas.gen.gz $B/gwas.sample \
        -M $B/genetic_map.txt --output-log $B/gwas.phased

#######################
## with a reference ###
shapeit -check -G $B/gwas.gen.gz $B/gwas.sample \
      --input-ref $B/reference.haplotypes.gz $B/reference.legend.gz $B/reference.sample \
        -M $B/genetic_map.txt --output-log $B/gwas.alignments

## extract out SNPs
awk '$1 == "Strand" {print $4}' $B/gwas.alignments.snp.strand | uniq -u > $B/SNPs_to_flip
awk '$1 == "Missing" {print $4}' $B/gwas.alignments.snp.strand | uniq -u > $B/SNPs_to_drop

more $B/SNPs_to_flip
more $B/SNPs_to_drop

############
## rerun ###
shapeit -check -G $B/gwas.gen.gz $B/gwas.sample \
      --input-ref $B/reference.haplotypes.gz $B/reference.legend.gz $B/reference.sample \
      
        -M $B/genetic_map.txt --output-log $B/gwas.alignments
