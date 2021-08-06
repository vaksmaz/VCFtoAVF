####	Variant pathogenicity determination

OVERVIEW: 

This is an in-silico tool that was designed to annotate variants and determine if they are potentially pathogenic.
To determine pathogenicity, variants are annotated with ANNOVAR and ClinVar, then ACMG criteria are applied using InterVar. 
InterVar results are adjusted by re-evaluating PVS1, PP5, the de novo status. We also adjust 
calls for genes which have special working groups recommendations including TP53, NF1 and a few others. 



EASY INSTALLATION:

For installation: 
Step 1:   
Install VEP somewhere

Step 2:
$: setup_vep_database.sh
### this will install everything in ~/.vep
OR 
Installing VEP database in other locations then ~/.vep:

$: cd /To/The/VEP/Database/
$: curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh37.tar.gz
$: curl -O http://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_refseq_vep_100_GRCh38.tar.gz

$: tar xzf homo_sapiens_refseq_vep_100_GRCh37.tar.gz
$:tar xzf homo_sapiens_refseq_vep_100_GRCh38.tar.gz

## if VEP is already installed do the following
$: mkdir ~/.vep
$: ln -S /path/to/vep/dataset/homo_sapiens_refseq ~/.vep/homo_sapiens_refseq

Step 3:
$: cd /to/VCFtoAVF/
$: setup.sh



#### To test the installation:
$: cd test/
$:  ./Test_run.sh





REQUIREMENTS 

1. PERL
2. Python 3.x (or Python2.7 if needed)
3. SNPEff 
4. VEP
5. ANNOVAR
6. InterVar

	3. SNPEff - Prior to running the VCFtoAVP please run the vcf through SNPEff ********



This code cluster houses ANNOVAR and InterVAR and those codes can be updated and put into the main folder.
For the most updated version of ANNOVAR go to (https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and
the most updated version of InterVar go to (https://github.com/WGLab/InterVar)


	4. VEP -- Install as required. 
	
If VEP is installed as a module, please load it prior to running. 
Install VEP databases as instructed than create a soft-link in the ~/.vep folder. 



TESTING AFTER INSTALLATION:

```

cd example
./Test_run.1.sh

```

Output:

1. filename.<date>.Vars_P_LP.avp    # All the variants classified as Pathogenic or Likely pathogenic 
2. filename.<date>.avp    # All the variants, regardless of status


RUNNING AVF:


```

cd example
/<PATH>/to/VCF_to_AVP_run.pl -vcf <input.vcf or vcf.gz> -r <HG19, GRCh37, 38, HG38 or GRCh38> or [-options ]

```



Running options:


    -vcf [file_name] - **Required  - vcf file - only required field. File must end in vcf, vcf.gz, txt or txt.gz
     
    -ids [id1,id2,id3] - All the IDs of the individuals as listed in the VCF. 
     		List in the order you want the data in. id2,id1,id4,id3
    
	-SubTot [integ] - total individuals in the study. default 1
     		if you want id 2 to be the proband and 1 to be parent list them as 
     		example  -SubTot 1 -ids id2,id1,id4,id3  - -- means 
     
    -rel - relation for each one, in the ID;s order
     	
    -b [string] - reference build, currently requires 19, HG19, GRCh37, 38, HG38, GRCh38. Default [hg19]
    
    -o [string] - output file name. default AVF_set.....
    
    -minDP [integer] - minimum depth for output.
    
    -DB [string] - Database for Annovar files. Default downloaded
    
    -Anno [string] - AnnoVar to be used. Default downloaded
    
    -tumor [string] - please indicate which is sample name is a tumor. 
    
    -AnnovFL [string] - use only if you already ran annovar and want to bypass it. This flag will stop Annovar from running.
    
    -IntervFL [string] - use if you already ran InterVar. This flag will stop intervar from running.
    
    -vep [string] - location of VEP
    
    -python [string] - python2.7 or python2.6










	
	
	
	
	
	
	
	
	
	
	
	
	
	