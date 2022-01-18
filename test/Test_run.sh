#!/bin/bash


../VCF_to_AVF_run_01.pl -vcf GMLOF.new.multiFam.snpeff.vcf \
	-b hg38 \
	-snpSep -snpSepFL family_file.txt \
	-short_fl \
	-o Test-families

##############################################

#### Testing a set without inputing any individuals
../VCF_to_AVF_run_01.pl -vcf GMLOF.new.multiFam.snpeff.vcf \
	-b hg38 \
	-snpSep -short_fl \
	-o Test-indv



##############################################

###  testing directly inputing indivudal IDs
../VCF_to_AVF_run_01.pl -vcf GMLOF.new.multiFam.snpeff.vcf \
	-b hg38 \
	-snpSep -short_fl \
	-ids Pro01,Father01,Mother01 \
	-rel Proband1,Father1,Mother1 \
	-o Single_fam_test_commLine_input






 