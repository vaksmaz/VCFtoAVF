#!/bin/bash

../VCF_to_AVF_run.pl -vcf GMLOF.new.vcf \
-ids BS_ER6ZP72Y,BS_PVCGJJ3D1,BS_C2GFCVS2 \
-rel proband,mother,tumor \
-tumor BS_C2GFCVS2 -o test.trio-zvtest \
-b hg38
