#!/usr/bin/perl


use POSIX qw(strftime);
my $today = strftime "%m-%d-%Y", localtime; # today's date-on output
use FindBin '$Bin'; # full path of bin -- $Bin


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR "
	This program is designed to convert individual vcfs to annotated file MAF equivelent format.
	We call this format AVF (annotatated Variant File). If you see this messege there may be an issue with the input files. 
	";
	
	print STDERR " Usage: $0  <your.pvs1> <your.avf.temp> \n\n";
	print STDERR "ERROR: One of the files required is missing - most likely not generated ";


	exit($return);
}


my $filename = $ARGV[1];  # autopvs1 file

if ($filename =~ /gz$/){
	open($fh, "gunzip -c $filename |") or die "gunzip $filename: $!";
}
else {
	open $fh, $filename || die "Can't read file head.txt file '$filename' [$!]\n";
}

while (<$fh>) {
	chomp;	
	
	@l = split "\t", $_;
	$l[0] = uc $l[0];
	$locTRG{$l[0]}{$ResPVS1} = $l[-1];
	$locTRG{$l[0]}{$criter} = $l[-2];
	
}

$filename = $ARGV[0]; # temp AVF file.


if ($filename =~ /gz$/){   open($fh, "gunzip -c $filename |") or die "gunzip $filename: $!";   }
else {   open $fh, $filename || die "Can't read file head.txt file '$filename' [$!]\n";   }

$filename =~ s/avf.tmp.gz$|avf.vcf.gz//;
$filename =~ s/avf.tmp$//;

$outfile = "$filename.$today.avf";
open(my $out, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";


$outfile1 = "$filename.$today.Vars_P_LP.avf";
open(my $out1, '>', $outfile1) or die "Can't read file '$outfile1' [$!]\n";


$header = <$fh>; chomp $header; $header =~ s/IntrvarRerun/IntrvarRerunRes	IntrvarRerunSupport/;	

print $out  "L-LP_Call_final	Reasoning_for_call	InterVar_recount_result	assertion_criteria	Recount_InterVar_support	PP5	PVS1_Correction	";
#print  $out "#VariationID	Hyperlink	Evidence_of_P/LP	Var_call_BadgeLab	Call_descrip	Alternant_flag_ExpertPanel	Alternant_flag_Badge	Alternant_flag_NonBadge	BadgeLabClinSig=Num	ReviewPanle=Num	NonBadgeLabClinSig=Num	$header\n";
print $out "$header\n";


print  $out1 "L-LP_Call_final	Reasoning_for_call	InterVar_recount_result	assertion_criteria	Recount_InterVar_support	PP5	PVS1_Correction	";
#print  $out1 "#VariationID	Hyperlink	Evidence_of_P/LP	Var_call_BadgeLab	Call_descrip	Alternant_flag_ExpertPanel	Alternant_flag_Badge	Alternant_flag_NonBadge	BadgeLabClinSig=Num	ReviewPanle=Num	NonBadgeLabClinSig=Num	$header\n";
print $out1 "$header\n";

while (<$fh>) {
	chomp;

	@l = split "\t", $_;	

	my @res;
	#my $temp = "NA";
	
	#my $clinvR = join "\t", @l[7..19];
	#$temp = "$_";


	my $return = rerun_InterVar($l[21],$_);
	my $ret1 = "$return	$_";
	
	my @res = split "\t", $ret1;
	my $calls = "NA";
	
	
	
	if ($res[15] =~ /Benign|benign|Uncertain|athogenic/) { ## check ClinVar
		$calls = "$res[15]	ClinVar $res[15] -- $res[16]";
	
	}
	
	elsif ($res[0] =~ /Benign|benign|Uncertain|athogenic/) { ## check InrterVar
		$calls = "$res[0]	InterVar_Recalculated $res[0]";
		
	}
	else {	
		$calls = "Uncertain_significance_noInfo	No clinVar or InterVar data";
	}
	
		
	if ($l[75] =~ /^splice|Splice/ and $res[15] !~ /Pathogenic|Likely pathogenic|Likely_pathogenic/) {		
		$calls = "VUS	***Splice_variant_without_clinVar_Support";
		#$calls = "$res[15]	***Splice_variant_without_clinVar_Support";
	}
	
	$calls =~ s/Likely pathogenic/Likely_pathogenic/;
	$calls =~ s/Uncertain significance/Uncertain_significance/;
	$calls =~ s/Likely benign/Likely_benign/;
	
	print $out "$calls	$ret1\n";
	#print "$calls	$ret1\n";

	if ($calls  =~ /^Pathogenic|^Likely pathogenic|^Likely_pathogenic/){
		
		print $out1 "$calls	$ret1\n";
		#print "$calls	$ret1\n";
	}
	
}

################################################################
##################################### SUBs
##############################################################


sub rerun_InterVar() {

	my ($INTV,$temp) = @_;
	$INTV =~ s/, /,/g;
	$INTV =~ s/^ //;
	@interv = split " ", $INTV;
	$PLP = "NA";
	$BLB = "NA";
	my %h;
	my $autoPV = "NA";
	my @ll;
	@ll = split "\t", $temp;
	my $ii;
	
	#print "ll[9]  --- $ll[9]\n";
	foreach $inf (@interv) {
		@inv = split '=', $inf;
		$inv[1] =~ s/\[|\]//g;
		@inv1 = split ',', $inv[1];
		$xi = 0;
		$num = 0;
		

		
		if ($inf =~ /PP/ and $ll[9] =~ /Evidence_of_P\/LP/) {
			$num = $inv1[4];
			$inv1[4] = 1;
			
			$ii = $inv1[4];
			
			

		}
		
		elsif  ($inf =~ /PP/ and $ll[9] !~ /Evidence_of_P\/LP/) {
			$num = $inv1[4];
			$inv1[4] = 0;
			$ii = $inv1[4];
			
		}
		
		
		foreach $tic (@inv1){
			$xi = $xi + $tic;
		}
		
		$h{$inv[0]}{$cnt} = $xi;

	}
	
	$PVS1 = $h{PVS1}{$cnt} ;
	$PS = $h{PS}{$cnt} ;
	$PM = $h{PM}{$cnt} ;
	$PP = $h{PP}{$cnt} ;
	$BP = $h{BP}{$cnt} ;	
	$BS = $h{BS}{$cnt} ;
	$BA1 = $h{BA1}{$cnt} ;
	
	
	$l[1] =~ s/chr//;
	$qrun = "$l[1]-$l[2]-$l[3]-$l[4]";
			
	
	#if (exists $locTRG{$qrun} and $locTRG{$qrun}{$ResPVS1} !~ /NF2|NF4|SS2|SS4|SS7|DEL3|DEL5|DEL9|DUP2|DUP4|DUP5|IC5/) {
	if ($locTRG{$qrun}{$criter} !~ /Unset/ and $locTRG{$qrun}{$ResPVS1} =~ /NF1|SS1|DEL1|DEL2|DUP1|IC1/) {
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Set";
		$PVS1 = 1;		
	}
	elsif ($locTRG{$qrun}{$criter} =~ /Unset/ and $locTRG{$qrun}{$ResPVS1} =~ /NF1|SS1|DEL1|DEL2|DUP1|IC1/) {
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
		#$autoPV = "Strong:Unset-PVS1";
		$PVS1 = 0;	
		$PS = $PS+1;		
	}
	elsif ( $locTRG{$qrun}{$ResPVS1} =~ /NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2/) {
		$PVS1 = 0;
		$PS = $PS+1;
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
	}
	elsif ($locTRG{$qrun}{$ResPVS1} =~ /NF6|SS6|SS9|DEL7|DEL11|IC3/) {
		$PVS1 = 0;
		$PM = $PM+1;
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
	}
	elsif ($locTRG{$qrun}{$ResPVS1} =~ /NF6|SS6|DEL7|DEL11|IC3/) {
		$PVS1 = 0;
		$PP = $PP+1;
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
	}	
	elsif ($locTRG{$qrun}{$ResPVS1} =~ /IC4/) {
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
		$PVS1 = 0;
		$PP = $PP+1;
	}
	elsif (exists $locTRG{$qrun} and $locTRG{$qrun}{$ResPVS1} !~ /NF2|NF4|SS2|SS4|SS7|DEL3|DEL5|DEL9|DUP2|DUP4|DUP5|IC5/) {
		$PVS1 = 0;
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
	}
	elsif (exists $locTRG{$qrun} and $locTRG{$qrun}{$ResPVS1} =~ /na|NF0/) {
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
		$PVS1 = 0;
		
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";

	}
	elsif (exists $locTRG{$qrun} ) {
		$autoPV = "$locTRG{$qrun}{$ResPVS1}:Unset-PVS1";
	}
	
	#else {$autoPV = "NA";}

	
	
=cut	
Pathogenic - Criteria 1
  (i) 1 Very strong (PVS1) AND 
  	(a) ???1 Strong (PS1???PS4) OR 
  	(b) ???2 Moderate (PM1???PM6) OR	 
 	(c) 1 Moderate (PM1???PM6) and 1 supporting (PP1???PP5) OR
 	(d) ???2 Supporting (PP1???PP5)	
=cut
	
	if ($PVS1 > 0 and ($PS > 0 or $PM > 1 or ($PM > 0 and $PP > 0) or $PP > 1)  ) { #criteria 1
		$PLP = "Pathogenic	criteria 1";
	
	}
	
=cut	
Pathogenic - Criteria 2
  (ii) ???2 Strong (PS1???PS4) OR

 Pathogenic - Criteria 3
  (iii) 1 Strong (PS1???PS4) AND
 	(a)???3 Moderate (PM1???PM6) OR
	(b)2 Moderate (PM1???PM6) AND ???2 Supporting (PP1???PP5) OR
	(c)1 Moderate (PM1???PM6) AND ???4 supporting (PP1???PP5)
=cut

	elsif ($PS > 1) {$PLP = "Pathogenic	criteria 2";}
	elsif ($PS > 0 and ($PM > 2 or ($PM > 1 and $PP > 1) or ($PM > 0 and $PP > 4)  )  )  {
		$PLP = "Pathogenic	criteria 3";	
	}
	


	
=cut	
Likely pathogenic 
	(i) 1 Very strong (PVS1) AND 1 moderate (PM1??? PM6) OR
	(ii) 1 Strong (PS1???PS4) AND 1???2 moderate (PM1???PM6) OR
	(iii) 1 Strong (PS1???PS4) AND ???2 supporting (PP1???PP5) OR
	(iv)  ???3 Moderate (PM1???PM6) OR
	(v) 2 Moderate (PM1???PM6) AND ???2 supporting (PP1???PP5) OR
	(vi) 1 Moderate (PM1???PM6) AND ???4 supporting (PP1???PP5)
=cut	
		
	elsif ($PVS1 > 0 and $PM > 0) {$PLP = "Likely_pathogenic	criteria 1"; }
	
	elsif ($PS > 0 and $PM > 0) {$PLP = "Likely_pathogenic	criteria 2"; }
	
	elsif ($PS > 0 and $PP > 1) {$PLP = "Likely_pathogenic	criteria 3"; }
	
	elsif ($PM > 2) {$PLP = "Likely_pathogenic	criteria 4"; }
	
	elsif ($PM > 1 and $PP > 1) {$PLP = "Likely_pathogenic	criteria 5"; }
	
	elsif ($PM > 0 and $PP > 3) {$PLP = "Likely_pathogenic	criteria 6"; }
	
=cut 
Benign
	(i) 1 Stand-alone (BA1) OR
	(ii) ???2 Strong (BS1???BS4)
=cut
	
	elsif ($BA1 > 0 or $BS > 1) {$BLB = "Benign	criteria 1 or 2"; }
	
=cut 
Likely Benign	
	(i) 1 Strong (BS1???BS4) and 1 supporting (BP1??? BP7) OR
	(ii) ???2 Supporting (BP1???BP7)
=cut

	elsif ($BS > 0 and $BP > 0) {$BLB = "Likely benign	criteria 1"; }
	elsif ($BP > 1) {$BLB = "Likely benign	criteria 2"; }
	
	else {$BLB = "Uncertain_significance	Criteria_unmet";}

=cut
Uncertain  significance
	(i) non of the criteria were met. 
	(ii) Benign and pathogenic are contradictory. 
=cut

	if ($PLP =~ /Likely_pathogenic|Pathogenic/ and $BLB =~ /NA/) {$outcome = $PLP;}
	elsif ($PLP !~ /Likely_pathogenic|Pathogenic/ and $BLB =~ /benign|Benign|Uncertain/) {$outcome = $BLB;}
	else {$outcome = "Uncertain_significance	Conflicing";}
	
	
	#print $out "$outcome	PVS1=$PVS1 PS=$PS PM=$PM PP=$PP BP=$BP BS=$BS BA1=$BA1	$inv1[4]	$_";
	$return = "$outcome	PVS1=$PVS1 PS=$PS PM=$PM PP=$PP BP=$BP BS=$BS BA1=$BA1	$ii	$autoPV";
	return $return ;
	
	
  }