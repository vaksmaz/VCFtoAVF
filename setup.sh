#!/bin/bash

mkdir humandb/

#### Get needed perl extensions

cpanm Getopt::Long
cpanm FindBin
cpanm LWP::Simple

./Run_CliVar_Re-annot.pl -build hg38 -get-URL
./Run_CliVar_Re-annot.pl -build hg19 -get-URL

git clone https://github.com/WGLab/InterVar.git
mv InterVar/* ./

wget https://omim.org/static/omim/data/mim2gene.txt
mv mim2gene.txt InterVardb/
mv mim2gene.txt intervardb/

#### download for annovarDB refGene,knownGene,ensGene

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar knownGene humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVer humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGeneWithVer humandb/


##### download from AnnoVar 1KG,exac03nontcga and Gnomad

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03nontcga humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03nontcga humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_genome humandb/


./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/

##### download from AnnoVar avsnp150, cosmic70 and dbnsfp35a

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35a humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp35a humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cosmic70 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic70 humandb/


#### for InterVar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar esp6500siv2_all humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp33a humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp33a humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/

#### download hg19 and hg38 from UCSC site

###############   Instruction: unzip the reference files and run 

#### please load Samtools 

## $: samtools faidx fasta hg38.fa  
###            and 
## $: samtools faidx fasta hg38.fa

wget -P pvs/data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -P pvs/data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

gunzip pvs/data/hg38.fa.gz
samtools faidx pvs/data/hg38.fa 

gunzip pvs/data/hg19.fa.gz
samtools faidx pvs/data/hg19.fa

## mkdir ~/.vep
## cd ~/.vep
## curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh37.tar.gz
## curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh38.tar.gz

## tar xzf homo_sapiens_refseq_vep_100_GRCh37.tar.gz
## tar xzf homo_sapiens_refseq_vep_100_GRCh38.tar.gz





