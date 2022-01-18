#!/usr/bin/perl


######################################
# Author: Zalman Vaksman
# Date: 2020
#
# Update ClinVar vcfs
######################################

#use strict;

use Cwd 'abs_path';
use FindBin '$Bin'; # full path of bin -- $Bin
use Getopt::Long qw(:config no_ignore_case);
#use LWP::Simple;
use Cwd;
use POSIX qw(strftime);
my $today = strftime "%m-%d-%Y", localtime; # today's date-on output

# Updating clivar criteria

############################### Set Conditions and running info

my $FullPath = abs_path($0); #path of file
my $bin = $Bin;
my $bt = "vcf_GRCh38";
my $Build = "hg38";
my ($helpFlag);
my $cmd;

GetOptions(
	"h|?|help"	=> \$helpFlag,

#	"f=i"	=> \$flankingLen,
	"build=s"	=> \$Build,    # genome build, 19, 37 or 38
	"cv=s"	=> \$ClinVcf,    #clinvar file
	"subSum=s"	=> \$subSum,    #clinvar file
	"get-URL"	=> \$geturl,    #submission summary file
	"GotFile"	=> \$ClinVarFile,    #submission summary file
	
) || help(1);

#if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }

######### get more vars
if ($Build =~ /38$/) { $bt = "vcf_GRCh38";}
elsif ($Build =~ /19$|37$/) { $bt = "vcf_GRCh37";}
my $ClinVcf = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/$bt/clinvar.vcf.gz";
my $subSum = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz";
######### get help

help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR "\n
######################################
# Author: Zalman Vaksman PhD
# Date: 2020
#
# Update ClinVar vcfs
######################################
	
	USAGE:
	
	To automate clinvar download for either HG19 or HG38 -- if build not specified HG38 is defualt
	
	/Path/To/Run_CliVar_Re-annot.pl -get-URL [optionl: < -b hg38 >]    ### builds recognized - HG19, HG38, hg19, hg38, 19, 38,


	__________________________________________
	To download specific files from the web please provide URLs as follows. -build is not required
	
	/Path/To/Run_CliVar_Re-annot.pl -get-URL < -b hg38 > -cv <URL ClinVar.vcf file> -s <URL submission summary file>


	__________________________________________
	If bypassing download, meaning the ClinVar and Submission summary files are available -- usage is as follows
	
	/Path/To/Run_CliVar_Re-annot.pl -GotFile -cv <Path/to/clinvar....vcf.gz> -s <Path/to/submission_summary.txt>
	
	
	";
	

	exit($return);
}

if ($geturl and $ClinVarFile) {
	$return = 0 if(!defined $return);
	print STDERR "CAN NOT use both get-URL and GotFile\n";
	print STDERR "-get-URL -- To download from web \n\t To manually add URL add -cv <file loc or http:loc/of/file>\n";
	print STDERR "\n Usage: $0 -get-URL -cv <http:clinvar file location>   : for specific location\n\n";
	print STDERR "\n Usage: $0 -get-URL    : the most recent version of ClinVar at $ClinVcf\n\n";
	
	exit($return);
}


if ($ClinVarFile and ($ClinVcf =~ /www|http/ or $subSum =~ /www|http/)) {
	$return = 0 if(!defined $return);

	print STDERR " \n
	Error  --- using the -GotFile flag you must provide th -cv clinvar file and the -s submission_summary file as below
	
	If bypassing download, meaning the ClinVar and Submission summary files are available -- usage is as follows
	/Path/To/Run_CliVar_Re-annot.pl -GotFile -cv <Path/to/clinvar....vcf.gz> -s <Path/to/submission_summary.txt>
	
	See /Path/To/Run_CliVar_Re-annot.pl -help
	";
	
	exit($return);
}


# Finding the required bed files
#if(!defined $bed_files){ print STDERR "\nNeed file info for the bed file files!!\n\n"; help(1); }
#my $fin_bed = "$bed_files".".bed";
#if(! -e $fin_bed){ print STDERR "\n[Error] Could not find the x.bed file '$fin_bed'.\n"; exit(1); }



############# run codes
my $bin = $Bin;

if ($Build =~ /38$/) {my $bt = "vcf_GRCh38";}
elsif ($Build =~ /19$|37$/) {my $bt = "vcf_GRCh37";}
my $cmd;

my $Sub = "$Bin/submission_summary.curr.txt.gz";
my $Sub1 = "$Bin/submission_summary";
my $clvar = "$Bin/ClinVar.$Build.curr.vcf.gz";

if ($geturl) {
	print "running now\n\n";

	$cmd = "wget $ClinVcf -O $clvar";
	print "$cmd\n";
	system("$cmd");
	

	$cmd = "wget $subSum -O $Sub";
	print "$cmd\n";
	system("$cmd");
	
}


# Reset_submission_summery_matrix.step1.pl - first part

$cmd = "$Bin/Reset_submission_summery_matrix.step1.pl $Sub1";
print "\n\n	Running command ...........\n	$cmd\n\n";
system("$cmd");


$cmd = "$Bin/ClinVarVCF_restructure_part2.pl /$Sub1.Re_assessed.$today.txt $clvar";
print "\n\n	Running command ...........\n	$cmd\n\n";
system("$cmd");

$cmd = "gzip $Sub.Re_assessed.$today.txt.ReCalled_final.txt";
#system("$cmd");

$cmd = "gzip $Sub.Re_assessed.$today.txt";
#system("$cmd");


print "###########################    Output   ################################


We completed the following functions:

1. You have the most recent ClinVar file 
2. You have the most recent Sumbission Summary file
3. The intermediate file for the Sumbission Summary 
4. The FINAL OUTPUT FILE (ClinVar File) --  
";












