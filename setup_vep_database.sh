#!/bin/bash



mkdir ~/.vep
cd ~/.vep
curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh37.tar.gz
curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh38.tar.gz

tar xzf homo_sapiens_refseq_vep_100_GRCh37.tar.gz
tar xzf homo_sapiens_refseq_vep_100_GRCh38.tar.gz


