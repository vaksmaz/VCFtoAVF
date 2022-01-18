#!/usr/bin/perl

use POSIX qw(strftime);
use FindBin '$Bin'; # full path of bin -- $Bin

my $today = strftime "%m-%d-%Y", localtime; # today's date-on output



########## compiling the Var_citation file into a matrix 
##########  Part 2 of the process - creating the hyperlinks and the citation matirix

##### required for running
		### 1.  var_citations.txt
		### 2. Output_step1  -- 



#my $filename = "submission_summary.txt.gz.Combined_called.12-07-2020.txt"; ## intervar file
my $filename = $ARGV[0]; ## var_citations.txt required for this part

if ($filename =~ /gz$/){
	open($fh, "gunzip -c $filename |") or die "Unable to gunzip $filename: $!";
}
else {
	open $fh, '<', $filename or die "Can't read file '$filename' [$!]\n";
}

print "Opening part2 var_citations file $filename ........\n";


my %h;

while (<$fh>) {
	chomp;
	if ($_ =~ /^#VariationID/) {
		$header = $_;
		next;
	}
	
	@l = split "\t", $_;
	$h{$l[0]}{$line} = $_;
}
	
	

#my $outfile = $ARGV[1];



my $file = $ARGV[1];


if ($file =~ /gz$/){
	open($fh, "gunzip -c $file |") or die "failed to gunzip $file  : $!";
}
else {
	open $fh, '<', $file or die "Can't read file '$file' [$!]\n";
}

$file =~ s/.curr.vcf.gz$//;
print "Opening ClinVar_VCF file $file ........\n";

my $outfile =  "$file.$today.Reassessed.txt";
open(my $out, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";

print "Opening part 2 output file $outfile ........\n";

while (<$fh>) {
	chomp;
	my $prev = "NA";
	if ($_ =~ /^#CHROM/) { $head = $_;	print $out "$head	Previous_CLNSIG_call	$header\n";	next; }
	if ($_ =~ /^#/) { 	print $out "$_\n"; next; }

	@l = split "\t", $_;
	
	@INFO = split ';', $l[7];
	foreach my $k (@INFO) {
		if ($k =~ /^CLNSIG/) {
			$k =~ s/^CLNSIG=//;
			$prev = $k;
		}
		
		if ($k =~ /^CLNSIG/) {
			$k =~ s/^CLNSIG=//;
			$prev = $k;
		}
	}
	
	if (exists $h{$l[2]}) {
			
		print $out "$_	$prev	$h{$l[2]}{$line}\n";
	}
	else {
		print $out "$_	$prev	Does Not Exists\n"; 
	#print  "$_	$prev	No Evidence Exists\n";
	}
	
	
}


