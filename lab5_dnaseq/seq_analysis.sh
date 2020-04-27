PS1='\W> '
## get VCF files
mkdir -p VCF

wget -P VCF/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget -P VCF/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

## get vcftools
## go here: https://sourceforge.net/projects/vcftools/files/latest/download
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure --prefix=/cloud/project/lab5_dnaseq/vcftools
make
make install
export PATH="$PATH:/cloud/project/lab5_dnaseq/vcftools/bin"

## need this too
export PERL5LIB=/cloud/project/lab5_dnaseq/vcftools/src/perl/

## go back out
cd /cloud/project
vcftools # works

## need tabix too
cd /cloud/project/lab5_dnaseq
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar vxjf htslib-1.9.tar.bz2
cd htslib-1.9    # and similarly for bcftools and htslib
./configure --prefix=$HOME
make
make install
export PATH="$PATH:/cloud/project/lab5_dnaseq/htslib-1.9"
tabix # works

## and samtools
cd /cloud/project/lab5_dnaseq
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar vxjf samtools-1.9.tar.bz2
cd samtools-1.9    # and similarly for bcftools and htslib
./configure --prefix=$HOME
make
make install
export PATH="$PATH:/cloud/project/lab5_dnaseq/samtools-1.9"
samtools

## and bcftools
cd /cloud/project/lab5_dnaseq
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9    # and similarly for bcftools and htslib
./configure --prefix=$HOME
make
make install
export PATH="$PATH:/cloud/project/lab5_dnaseq/bcftools-1.9"
bcftools

## get SRA toolkit
cd /cloud/project/lab5_dnaseq
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar vxzf sratoolkit.current-centos_linux64.tar.gz
cd sratoolkit.2.10.5-centos_linux64
export PATH="$PATH:/cloud/project/lab5_dnaseq/sratoolkit.2.10.5-centos_linux64/bin"

## new in 2.10
vdb-config -i 

## test
fastq-dump -X 5 -Z SRR390728
## works


###############################
## examples using 1000G data ##
###############################
cd /cloud/project/lab5_dnaseq


## motivated by: https://vcftools.github.io/man_latest.html
## http://samtools.github.io/hts-specs/VCFv4.2.pdf

vcf=/cloud/project/lab5_dnaseq/VCF/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

## frequeny
vcftools --gzvcf $vcf --freq --chr 21 --out chr21_freq
## hmm this could take a while

## maf filter 
vcftools --gzvcf $vcf  --maf 0.05 --remove-indels --recode --recode-INFO-all --out /cloud/project/lab5_dnaseq/VCF/chr21_common
gzip VCF/chr21_common.recode.vcf

## r2 of genotypes
vcftools --gzvcf VCF/chr21_common.recode.vcf.gz  --geno-r2  --out VCF/chr21_common

## make small
tabix -h $vcf 21:9411239-11411239 > VCF/chr21_small.vcf
vcftools --vcf VCF/chr21_small.vcf --freq --chr 21 --out chr21_small_freq

## or using vcftools
vcftools --gzvcf $vcf --chr 21 --from-bp 9411239 --to-bp 11411239 --recode --recode-INFO-all > VCF/chr21_small_2.vcf

## depth
vcftools --vcf VCF/chr21_small.vcf --depth -c > depth_summary.txt

###################
## file examples ##

############################
## get reads from 1000G ####
# http://www.internationalgenome.org/category/fastq/
# http://www.internationalgenome.org/faq/where-are-your-sequence-files-located/

## this sample: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096
## http://www.internationalgenome.org/data-portal/sample/HG00096

mkdir -p /cloud/project/lab5_dnaseq/FASTQ

## just write out the first 100000 reads
fastq-dump -X 100000 -O /cloud/project/lab5_dnaseq/FASTQ/ --split-3 SRR062634
wc -l FASTQ/SRR062634_1.fastq
head FASTQ/SRR062634_1.fastq

#### BAMs ####
mkdir -p BAM
wget -P BAM/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam*

## view header
samtools view -H BAM/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

## view alignments
samtools view BAM/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam | head

## view stats
more BAM/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bas

## call variants, basic
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz
sed -i '1c\>20' chr20.fa 
samtools faidx chr20.fa

samtools mpileup -uf chr20.fa -AB -q0 -r 20:1000000-1100000 BAM/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam -o tmp.vcf
bcftools call -mv -O tmp.vcf > VCF/HG00096_chrom20_subregion.vcf.gz
## workflow examples ##
## http://www.htslib.org/workflow/#mapping_to_variant