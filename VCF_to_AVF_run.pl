#!/usr/bin/env perl

use Getopt::Long qw(:config no_ignore_case);
use Cwd;
use FindBin '$Bin'; # full path of bin -- $Bin


my $cwd = getcwd();

############################### Set Default Conditions -- and running info

my $name = "NA";
my $build = "hg38";
my $subNum = 1;
my $rel = "Proband";
my $output = "AVF_set";
my $build = "hg19";
my $annovar = "$Bin/table_annovar.pl";
my $intervar = "$Bin/Intervar.py";
my $database = "$Bin/humandb/";
my $vcfAnnov = "NA";
my $vcfInterV = "NA";
my $vep = "vep";
my $python = "python";
#my $python2 = "python";
my $dbnsfp = "dbnsfp35a";

GetOptions(
	"h|?|help"	=> \$helpFlag,
	"vcf=s"	=> \$vcf, # Normal VCF
	"ids=s" => \$name, #ids of all the individuals tested, comma separated
	"b=s"	=> \$build, 
	"SubTot=i" => \$subNum, # total numb of indi, default = 1
	"rel=s" => \$rel, # relationship
	"o=s" => \$output, #output
	"minDP=i" => \$minDepth,
	"DB=s" => \$database, # annovarDB
	"tumor=s" => \$tumor, # annovarDB
	"AnnovFL=s" => \$vcfAnnov, # annovarDB
	"IntervFL=s" => \$vcfInterV, # annovarDB
	"clinvarL=s" => \$clinv, # annovarDB
	"vep=s" => \$vep, # annovarDB
	"python=s" => \$python, # python2
#	"python2=s" => \$python2, # python2
	"dbnsfp=s" => \$dbnsfp, #dbnsfp version
) || help(1);



sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR "This program is designed to convert individual vcfs to annotated file MAF equivelent format.
	We call this format AVF (annotatated Variant File). This program annotates everything using annovar 
	and is currently set up to use HG19 or HG38. A more formal tool to use with other annotation references 
	is currently in progress. Also, this has to be used on respublica, again a more formal method will be completed in time.
	
	This software requires other programs - ANNOVAT and InterVar
	
	To update ANNOVAR and InterVar go to their respective websites. 

	";

	print STDERR " Usage: $0  -vcf <vcf>  options [ -m <min depth> -s <Normal/Tumor>]ll command> ]\n\n";
	print STDERR "  ex 1) $0 -vcf your.vcf -tumor your.tumor.vcf -b HG38";
	#print STDERR "  ex 2) $0 -o NBL_cases -b All_NBL_Omni -p /home/dir/plink1.9/plink -hm /5home/dir/hapmap_exp/hapmap.hg19.consensus.qc -a /home/dir/hapmap_exp/hapmap_ancestry_info.txt\n";


	print STDERR "
		------------------------------------
		Options  \n
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
    
    -python [string] - python3 
    
   ***** -python2 [string] - python2.7 or python2.6
    
    -dbnsfp [string] - version of ANNOVAR dbnsfp. Default dbnsfp35a

				\n\n";


	exit($return);
}

##################### put in conditions for error.

if (!defined $vcf){ print STDERR "\nNeed file info for the vcf file files!!\n\n"; help(1); }

#if(! -e $vcf){ print STDERR "\n[Error] Could not find the x.bed file '$vcf'.\n"; exit(1); }



####################### split multi-non reference allele info





#### make tem directory
my (@l, $fold);

	$fold = $vcf;
	$fold =~ s/.txt.gz$|.vcf.gz$|.vcf$|.txt$/_TMPfolder/;
	
	system("mkdir $fold/");
	system("cp $vcf $fold/");
	chdir "$fold";

	
##################### Set individuals Test

if ($vcf =~ /gz$/){
	open($fh1, "gunzip -c $vcf |") or die "gunzip $vcf: $!";

}
elsif ($vcf =~ /vcf|txt/) {
	open $fh1, $vcf || die "Can't read file head.txt file '$vcf' [$!]\n";

}
else {print "ERROR: need a VCF or TXT file -- please make sure the file is ended in vcf, txt, vcf.gz or txt.gz";
		help(1);
	}


my $outfile_INTV = $output . ".vcf.InterVar.tmp";
open(my $out2, '>', $outfile_INTV ) or die "Can't read file '$outfile' [$!]\n";


my $outfile = $output . ".vcf.tmp";
open(my $out1, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";

my ($vcf_info, $snpeff_info, @altAlltmp);


while (<$fh1>) {
	chomp;
	@l = split "\t", $_;
	if ($_=~ /INFO=<ID=ANN/) {chomp; print $out1 "$_\n";
		$_ =~ s/\|/\t/g; $_ =~ s/ //g; 
		@l = split '\'', $_;
		$snpeff_info = $l[1];
	}	 #get snpeff header info from vcf

	elsif ($_=~ m/#CHROM/) {  
			$vcf_info = $_;
			print $out1 "$_\n";
			print $out2 join "\t", @l[0..7], "FORMAT\tGeno\n";
		} #vcf header
	
	elsif ($l[4] =~ /,/) {
		@altAlltmp = split ",", $l[4];
		$l[7] = "MULTIALL=$l[4];$l[7]";
		
		foreach my $Aall (@altAlltmp) {
		
			if ($Aall eq "*") {
				next;
			}
			else {
				$l[4] = $Aall;
				print $out1 join "\t", @l, "\n";
				print $out2 join "\t", @l[0..7], "GT\t0/1\n";
				}
			}			
		}
	elsif ($_ =~ /##/) {
		print $out1 "$_\n";
		print $out2 "$_\n";
	
	}
		else {
			print $out1 "$_\n";
			print $out2 join "\t", @l[0..7], "GT\t0/1\n";
			}
		
}

close $fh1;
close $out1;
close $out2;
	####### Snpeff annotation error code
if (!defined $snpeff_info ) {
	print "\n\n********  PLEASE RUN AN ANNOTATION CODE SUCH AS SnpEff or VEP:\nhttps://pcingola.github.io/SnpEff/  ********\n\n";
	help(1);
	print "\n\n ---Error: No SnpEff annotation: \n###### PLEASE RUN AN ANNOTATION CODE SUCH AS SnpEff or VEP:\nhttps://pcingola.github.io/SnpEff/ -----  \n\n";
	die;	
	}





##################### Set individuals Test
my $c = 0;
my @vinf;
my $cases = "Case_notes";
if ($name eq "NA") { 
	@vinf = split "\t", $vcf_info;
	$name = $vinf[-1];}

my @IDs = split ",", $name;
my @relat = split ",", $rel;
my (@IDset, %tmrHS);
my @tmr = split ',', $tumor; 
foreach (@tmr) { $tmrHS{$_} = $_; }

foreach (@IDs) {
	if (exists $tmrHS{$_} ) {
		$cases = $cases."	Tumor_ID	Tumor_Notes	Tumor_Allele1	Tumor_Allele2	genotype	Tumor.Allele_1_Depth	Tumor.Allele_2_Depth";	
	}
	else {
		$c++;
		$cases = $cases."	Sub$c.ID	Sub$c.note	Ind.$c.Normal_Allele1	Ind.$c.Normal_Allele2	Ind.$c.genotype	Ind.$c.Allele_1_Depth	Ind.$c.Allele_2_Depth";
	
	}	
		
}


###################### Run Annovar


if ($build =~ /19/ or $build =~ /37/ ) { $bld = "hg19";}
else { $bld = "hg38";}



if ($vcfAnnov eq "NA") {
	my $cmd = "perl $annovar $outfile $database -buildver $bld "
	. "-protocol refGene,knownGene,ensGene,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas"
	. ",1000g2015aug_afr,1000g2015aug_all,exac03nontcga,gnomad211_genome"
	. ",avsnp150,cosmic70,$dbnsfp"
	. " -operation g,g,g,f,f,f,f,f,f,f,f,f,f --out $output -nastring . -vcfinput";
	print "Running Annovar ------------ \n\n$cmd\n\n";
	system($cmd);
	
	$vcfAnnov = "$output.$bld"."_multianno.txt";
}	




print "running cleanup ................\n\n";
my $cmd = "rm $output*dropped";
system("$cmd");

my $cmd = "rm $output*filtered";
system("$cmd");

my $cmd = "rm $output*log";
system("$cmd");

my $cmd = "rm $output*function";
system("$cmd");


######################################## run InterVar

if ($vcfInterV eq "NA") {
	$cmd = "$python $Bin//Intervar.py -b $bld -i $outfile_INTV "
	. " --input_type=VCF -o $output.intervFL --table_annovar=$Bin/table_annovar.pl "
	. "--convert2annovar=$Bin/convert2annovar.pl "
	. "--annotate_variation=$Bin/annotate_variation.pl "
	. "-d $Bin/humandb "
	. "-t $Bin/intervardb";
	
	print "$cmd\n\n";
	system($cmd);
	
	$vcfInterV = "$output.intervFL.$bld"."_multianno.txt.intervar";	
}



##################################  process InterVar File


open $fh, $vcfInterV || die "Can't read file head.txt1 file '$vcfInterV' [$!]\n";

print "\nOpening intervar data $vcfInterV\n\n";
my(%intv);

while (<$fh>) {
	chomp;
	@l = split "\t", $_;
	$itv = "chr$l[0]-$l[1]-$l[3]>$l[4]";
	$l[13] =~ s/ InterVar: //;
	$l[13] =~ s/ PVS1/\tPVS1/;
	$intv{$itv}{$call} = $l[13]; # putting result into a memory

}
close $fh;



########################### Processing ClinVar

if ($clinv) {next; }
else {
	opendir my $dh, $Bin or die "Could not open '$Bin' for reading '$!'\n"; # open dir
	
	while (my $thing = readdir $dh) {
		if ($thing =~ /ClinVar/  and $thing =~ /Reassessed.txt/ and $thing =~ /$bld/) {
        	$clinv = "$Bin/$thing";  ## same as #$clinv = `ls $Bin/ClinVar.$bld.*.Reassessed.txt*`; print "$clinv\n";}
    	}	 	
	}
}

close $dh;

print "ClinVar file used -- $clinv\n";

if ($clinv=~ /gz$/){
	open($fh, "gunzip -c $clinv |") or die "gunzip $clinv: $!";
}
else {
	open $fh, $clinv || die "Can't read file clinvar file '$clinv' [$!]\n";
}

#my $testset = "test.output";
#open(my $tst1, '>', $testset) or die "Can't read file 'test.output' [$!]\n";


$header = <$fh>; chomp $header;

my ($clinLoc,$clinCall,$calls);
my %clinHash;
while (<$fh>) {
	chomp;
	#print "$_\n";
	$_ =~ s/^chr//;
	@l = split "\t", $_;
	$clinLoc = "chr$l[0]-$l[1]-$l[3]-$l[4]";
	$clinCall = join "\t", @l[9..21];
	$clinHash{$clinLoc}{$calls} = $clinCall;
	
}

close $fh;



################## VarMod DBset

my (%GeneBasDB, %selecNM, %standNM);
my $dbVarMod  = "$Bin/VarMod/Var_info.db";
open $fh, $dbVarMod || die "Can't read file head.txt3 file '$filename' [$!]\n";

while (<$fh>) {
	chomp;
	@l = split "\t", $_;
	
	if ($_ =~ /Validated/ and $l[4] =~ /NM/ and $l[5] =~ /NP/) {
		$standNM{$l[2]}{NP} = $l[5];
		$standNM{$l[2]}{NM} = $l[4];
	
	}
	elsif ($l[4] =~ /NM/ and $l[5] =~ /NP/ and $l[1] =~ /Reviewed/) {
		$selecNM{$l[2]}{NM} = $l[4];
		$selecNM{$l[2]}{NP} = $l[5];
	}
	
}



########################## Creating the Annotated File

$outfile = "$output.$bld.avf.tmp";
open($out, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";


my $outfilePVS = $vcf . ".tmp";
open(my $outPVS, '>', $outfilePVS) or die "Can't read file '$outfile' [$!]\n";

print $outPVS "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	subjInfo\n";



open $fh1, $vcfAnnov || die "Can't read file head.txt4 file '$vcfAnnov' [$!]\n";
print "opening $vcfAnnov \nOpening $outfile ......................\n\n";


my $chng = "Otherinfo1	Otherinfo2	Otherinfo3	$vcf_info";
my $header = <$fh1>;
$header =~ s/Otherinfo/$chng/;
$header =~ s/\./_/g;
$header =~ s/\n//g;

	#print "$header\n";
###################### print header

# print ClinVar header Info
my $clinvHeader = "#VariationID	Hyperlink	Evidence_of_P/LP	New_ClinSig_Call	Call_descrip	Alt_ClinSig_call	Alt_Call_descrip	Alternant_flag_ExpertPanel	Alternant_flag_Badge	Alternant_flag_NonBadge	BadgeLabClinSig=Num	ReviewPanle=Num	NonBadgeLabClinSig=Num";
#my $clinvHeader = "#VariationID	Hyperlink	Evidence_of_P/LP	New_ClinSig_Call	Call_descrip	Alt_ClinSig_call	Alt_Call_descrip	Alternant_flag_ExpertPanel	Alternant_flag_Badge	Alternant_flag_NonBadge	BadgeLabClinSig=Num	ReviewPanle=Num	NonBadgeLabClinSig=Num";

# print InterVar and the rest of the header Info

print $out "Hugo_symbol	Chromosome	Loc	Reference	Alternant_alleles	Overall_Variant_Type	LocusAlleles	$clinvHeader	"
		.	"InterVar_res	InterVar_Support	$cases	"
		.	"De_novo	De_novo_info	LOF_Y_N	LOF_info	NMD_Y_N	NMD_info	$snpeff_info	$header\n";
		
	$snpeff_info =~ s/\./_/g;	
#print $outPVS "#Uploaded_variation	SYMBOL	Feature	CANONICAL	PICK	Consequence	HGVSc	HGVSp	HGVSg	EXON	INTRON\n";
#######################

#print "$header\n";

my ($xyz,$xyy) = 0; # Set counters


############################# creating AVF file


while (<$fh1>) {
		chomp;

	
	##################### counting how many snps have been converted
	$xyz++;
	if ($xyz == 10000) { $xyy++; $totnum_count = $xyz*$xyy;
		print "completed $totnum_count\n"; $xyz=0;
		}


	############################# Zeroing out all necessary hashes and arrays
	my (%info, @eff, @snp, %h, %dept, @alleles, @c, @tmp1, %snpE );
    my ($idx, $idy) = 0; # new counters
    my ($Variant_Type, $id, $all_alleles, $tmpINF, $itvs, $IDclinv, $currentALT);

	$_ =~ s/\|\,/\|NA\,/g;
	$_ =~ s/\|\t/\|NA\t/g;
	@snp = split '	', $_;
	

	######### set all the hashes needed in one shot.

	%h = map { $_ => $snp[$idx++]} split ( "\t", $header ); # map snp fields into a hash by header
	my @c = split /\;/, $h{INFO};

	####### Correcting fields in the snpEff results  ###########

    foreach $field (@c) {	if ($field !~ /=/) {$field  = $field . "=.";} 	else {$field;} }
	$h{INFO} = join "\;", @c;

	
    $tmpINF=$h{INFO};
	#	%info = map{split /=/, $_ }(split /\;|DB\;|POSITIVE_TRAIN_SITE\;/,$tmp ); # hash the info, splitting by ';'
	%info = map{split /=/, $_ }(split /;/,$tmpINF ); # hash the info, splitting by ';'
	$currentALT = $h{ALT};
	if (exists $info{MULTIALL} ) {$h{ALT} = $info{MULTIALL};}
	
#print "$info{DP}	$info{MULTIALL}	$h{ALT}		currentAlt-- $currentALT\n";
	
#print "$currentALT	$info{MULTIALL}\n";   ####################
	#  $info{MULTIALL}
	
	my @eff1 = split ',', $info{ANN}; 
	my $ind_snp;
	my @eff_q;
	foreach $ind_snp (@eff1) {
		$ind_snp =~ s/$/-INFO/;
		@eff_q = split "\|", $ind_snp;

		if ($eff_q[0] eq $currentALT) {
			@eff = @eff_q;
			last;
		}	
	}
	
	
	
	
	
	$eff1[0] =~ s/$/-INFO/; my @eff = split '\|', $eff1[0]; # split into the snpEff results by '|'
	foreach my $str (@eff) { if ($str eq '') {$str = "."; }  }	
	my $idy = 0;
	%snpE = map { $_ => $eff[$idy++]} split ( "\t", $snpeff_info );
	 
	
	
	######################################################################################
	
	
	
	######################   print Gene, loc, ref and alt alleles

	my $it = "$h{Chr}-$h{POS}-$h{Ref}>$h{Alt}";  #info for InterVar	
	print $out "$h{Gene_refGene}	$h{Chr}	$h{POS}	$h{REF}	$currentALT	"; # print Gene, loc, ref and alt alleles	
	#print $out "$h{Gene_refGene}	$h{Chr}	$h{POS}	$h{REF}	$h{$ALT}	"; # print Gene, loc, ref and alt alleles	
	#######################		
	
	if ($h{Ref} =~ /-/){ $Variant_Type = "INS";}
	elsif ($h{Alt} =~ /-/){ $Variant_Type = "DEL";}
	elsif ($h{Alt} =~ /0/ or $h{Alt} eq '*') {$Variant_Type = "Remove - middle of deletion";}
	else {$Variant_Type = "SNP"; }
	print $out "$Variant_Type	$it	";  # print variant type (snp, INDEL)
		
	###################### Add ClinVar INFO
	$h{Chr} =~ s/^chr//;
	#$IDclinv = "chr$h{Chr}-$h{POS}-$h{REF}-$h{ALT}";
	$IDclinv = "chr$h{Chr}-$h{POS}-$h{REF}-$currentALT";

	if (exists $clinHash{$IDclinv}) {
		### $clinHash{$clinLoc}{$calls}		
		print $out "$clinHash{$IDclinv}{$calls}\t";	
#print $tst1 "New output 		$IDclinv\n";
	}
	else {
		print $out $temp = "NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	";	
	}

	
		
	######################   Add InterVar - ran here not DB
	$h{Chr} =~ s/chr//;
	
	$itvs = "chr$h{Chr}-$h{Start}-$h{Ref}>$h{Alt}";
	$itvs1 = "chr$h{Chr}-$h{Start}-$h{Ref}>0";
	#print "$itvs\n";
	if (exists $intv{$itvs}  ) {
		print $out "$intv{$itvs}{$call}	$rel	";
	}
	elsif(exists $intv{$itvs1}) {
		print $out "$intv{$itvs1}{$call}	$rel	";
	}
	else {print $out ".	.	$rel	";}
	
	
	######################   # print variant type (snp, INDEL)

	######################
	
	
	############################### breakdown of all alleles
	$c = 0;
	$all_alleles = "$h{REF},$h{ALT}";
	#print "ALTallele-alleles $h{ALT}		All alleles -- $all_alleles\n";	
	#$all_alleles = "$h{REF},$h{ALT}";
	@alleles = split ',', $all_alleles; # all the alleles in for this locus.

	@qqq = split "\t", $header;

	foreach $id (@IDs) {
		$idy=0;
		
		@allelect = split '\:', $h{$id}; #Set sample info parent2 -- allele numb	
		%dept = map { $_ => $allelect[$idy++]} split( '\:', $h{FORMAT} ); # map sample info to hash for the indiv
		
		@tmp = split '/', $dept{GT}; # get genotype info
		@tmp1 = split ',', $dept{AD}; # get allele depth at each allele
		#print "ID = $id		$h{FORMAT}	$h{$id}\n";
		############################  Print individual and allele info for each individual
		print $out "$id	$relat[$c]	$alleles[$tmp[0]]	$alleles[$tmp[1]]	$dept{GT}	$tmp1[$tmp[0]]	$tmp1[$tmp[1]]	";
		#############################
		$c++;
		
		
	}
	
	#################### De novo - find if the var is de novo in trio data
	
	if (exists $info{loConfDeNovo}) {print $out "LowConfidence	$info{loConfDeNovo}	";}
	elsif 	(exists $info{hiConfDeNovo}) {print $out "HighConfidence	$info{hiConfDeNovo}	";}
	else {print $out "None	.	";}
	
	##################### 
		

	# snpEff_LOF_Y_N
	if (exists $info{LOF}) {print $out "Y\t$info{LOF}\t";}
	else {print $out "N\t.\t";}
		
	
	# snpEff_NMD_Y_N # nonsense mediated decay
	if (exists $info{NMD}) {print $out "Y\t$info{NMD}\t";}
	else {print $out "N\t.\t";}
	
	
	
	print $out join ("\t", @eff);
	my @headr = split "\t", $header;
	my $sampInfo;
	foreach $sampInfo (@headr) {print $out "\t$h{$sampInfo}";} 
	#print $out "$_";	# print multianno file info -

	print $out "\n";	#end of line
	
	
	if ($h{INFO} =~ /LOF/ or $h{INFO} =~ /NMD/ or $intv{$itvs}{$call} =~ /PVS1=1/ or $h{INFO} =~ /splice|Splice|stop|Stop|start|frameshift|HIGH|MODERATE/) {
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	subjInfo
		if ($h{Chr} =~ /chr/) {
			print $outPVS "$h{Chr}	$h{POS}	.	$h{REF}	$currentALT	$h{QUAL}	$h{FILTER}	DP=stand	GT	0/1\n";
		}
		else {
			print $outPVS "chr$h{Chr}	$h{POS}	.	$h{REF}	$currentALT	$h{QUAL}	$h{FILTER}	DP=stand	GT	0/1\n";
		}
		
	}

	
}


close $out;
close $fh1;
close $autoPVS;

%clinHash = undef;
%intv = undef;

my $refbld;
my $grch;
if ($build =~ /19/ or $build =~ /37/ ) { $refbld = 19; $grch = 37;}
else { $refbld = 38; $grch = 38;}




$cmd = "$vep --offline --refseq --use_given_ref " .
    "--assembly \"GRCh"."$grch\" " .
    "--fork 4 " .
    "--canonical " . 
    "--flag_pick " .
    "--hgvs --hgvsg --symbol " .
    "--distance 500 " .
    "--exclude_predicted " .
    "--numbers " .
    "--lookup_ref " .
    "--input_file $vcf.tmp " .
    "--output_file $vcf.vep --force_overwrite --no_stats " .
    "--fasta $Bin/pvs/data/hg$refbld.fa " .
    "--tab --fields \"Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,Consequence,HGVSc,HGVSp,HGVSg,EXON,INTRON\" ";
    
print "Running VEP\n$cmd\n\n";
system("$cmd");




$cmd = "$Bin/grab_Pick1.pl $vcf.vep > $vcf.1.vep";
print "Running VEP\n\n";
system("$cmd");


$outfile = "$Bin/pvs/config.ini";
open($out, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";



print $out "
[DEFAULT]

ref = $Bin/pvs/data/hg$refbld".".fa
trans = $Bin/pvs/data/refGenePlus_20191020.gpe
domain = $Bin/pvs/data/PM1.domain.bed
hotspot = $Bin/pvs/data/PM1.hotspot.bed
curated_region = $Bin/pvs/data/PM1.expert_curated.bed
pathogenic_ref = $Bin/pvs/data/clinvar_pathogenic_20200106.vcf
pvs1levels = $Bin/pvs/data/PVS1.level
exon_lof_popmax = $Bin/pvs/data/exon_lof_popmax.bed
gene_trans = $Bin/pvs/data/clinvar_trans_stats_20200106.tsv
gene_alias = $Bin/pvs/data/hgnc.symbol.previous.tsv

";
close $out;



$cmd = "$python $Bin/pvs/vep_lof_filter.py $vcf".".1.vep";
print "Running calculating ...........\n\n";
system("$cmd");


$cmd = "$python $Bin/pvs/Run_pvs_res.1.py $vcf".".1.vep.lof $Bin > $vcf".".1.vep.pvs1";
print "Running Final results for pvs correction\n\n";
system("$cmd");



$cmd = "$Bin/InterVar_recalculate.pl $output.$bld.avf.tmp $vcf".".1.vep.pvs1";
print "\n\n$cmd\n\n";
system("$cmd");



##########################################################################

###############  Clean set



print "Cleaning up the temp files\n\n";


$cmd = "mv *avp ../";
system("$cmd");




system("rm *filtered");
system("rm *dropped");
system("rm *function");
system("rm *log");

=cut



chdir ($cwd);

$cmd = "rm -rf $fold";
system("$cmd");


print "Cleaning up the temp files\n";
system("rm *filtered");
system("rm *dropped");
system("rm *function");
system("rm *log");
system("rm *avinput");
system("rm *tmp");
system("rm *vcf");
system("rm *vcf.gz");
system("rm *vep");
#system("rm *all");
#system("rm *lof");
#system("rm *pvs1");
system("rm *fa");
#system("gzip *");












