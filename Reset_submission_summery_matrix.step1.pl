#!/usr/bin/perl

############################## code changes and issues and their dates. 
##########  11-24-2020 -- Discussion of variants with CBTTC found issues with calls that have no Badge labs but nonbadge lab info. -- Rule had been changed 
##########  12-3-2020 -- Previous rule changed, to now flag those and create alternant call column.  Rule had been changed 
##########   12-3-2020 -- As per Sharon, Jung and Doug's request 12-3-2020 - We keep the original call for rule 2 but flag and provide an ALT call. 
##########   12-4-2020 -- 12-3-2020 call changes are now reverted back and have an alt call column and a flag. 



use List::Util qw(sum);
use POSIX qw(strftime);
use FindBin '$Bin'; # full path of bin -- $Bin

my $today = strftime "%m-%d-%Y", localtime; # today's date-on output

%month2numb = qw(
    Jan 01  Feb 02  Mar 03  Apr 04  May 05  Jun 06
    Jul 07  Aug 08  Sep 09  Oct 10 Nov 11 Dec 12
);


############ Collect ClinVar submission_summary _new information
##	252820	Sep 16, 2011:-:Breast-ovarian cancer, familial 2:C2675520:Breast-ovarian cancer, familial 2:no assertion criteria provided:clinical testing:germline:1:Sharing Clinical Reports Project (SCRP):SCV000297423.1:BRCA2;
##	210186	Jun 27, 2014:-:Spinocerebellar ataxia, autosomal recessive 10:C3150998:Spinocerebellar ataxia, autosomal recessive 10:criteria provided, single submitter:clinical testing:germline:na:Genetic Services Laboratory, University of Chicago:SCV000246417.1:ANO10;

#my $filename = "/Users/vaksmanz/Desktop/database_info/NewClinVar/submission_summary.txt.gz"; ## intervar file




$filename = "$Bin/VarMod/badgeList";  #Set Badge labs

$badgeList = split "\n", `cat $Bin/VarMod/badgeList`;
$badgeList =~ s/\n//g;




my $filename = $ARGV[0]; #Submission Summary file 

my $outfile =  "$Bin/submission_summary.Re_assessed.$today.txt";
open(my $out, '>', $outfile) or die "Can't read file '$outfile' [$!]\n";



if ($filename =~ /gz$/){
	open($fh, "gunzip -c $filename |") or die "gunzip $filename: $!";
}
else {
	open $fh, $filename || die "Can't read file '$filename' [$!]\n";
}

print "Opening submission_summary file $filename ........\n";
my %h;

my (@l,@info,@date);
my ($month,$pdate,$cl);

	##  Collect data -
	## data obtained - date|ClinSig|Submitter|ReviewStatus|ReportedPheno|SCV|Gene|Description



print "Opening output file $outfile ........\n";

print $out "#VariationID	Hyperlink	Evidence_of_P/LP	New_ClinSig_Call	Call_descrip	Alt_ClinSig_call	Alt_Call_descrip	Alternant_flag_ExpertPanel	Alternant_flag_Badge	Alternant_flag_NonBadge	BadgeLabClinSig=Num	ReviewPanle=Num	NonBadgeLabClinSig=Num	date|ClinSig|Submitter|ReviewStatus|ReportedPheno|SCV|Gene|Description\n";

my $LineCounts = 0;
while (<$fh>) {
	chomp;	
	$LineCounts++;
	if ($LineCounts == 10000) {
		print ".";
		$LineCounts=0;
	}
	
	if ($_ =~ /^#/){
		next;
		}
	
	@l = split "\t", $_;
	
	
			
		@info = split ':', $l[2];
		@date = split /,| /, $info[0];
		$month = $month2numb{$date[0]};
		$pdate = "$date[3]"."$month"."$date[1]";  ## collect date
	## data obtained - date|ClinSig|Submitter|ReviewStatus|ReportedPheno|SCV|Gene|Description
	
		$cl = "$pdate|$l[1]|$l[9]|$l[6]|$l[5]|$l[10]|$l[11]|$l[3]"; ## ClinVar _info

		$h{$l[0]}{$cl} = $cl;	

}

close $fh;

foreach $k ( keys %h ) {
		my @var_info;
		my $i = 0;
		my $badge_info = "NA";
		print $out "$k	=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/clinvar/variation/$k/\",\"NCBI_$k\")	";

		#print "$k	=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/clinvar/variation/$k/\",\"NCBI_$k\")	\n";

		
	foreach $m ( keys %{$h{$k}} )  {
		$i++;
		push @var_info, $h{$k}{$m};
	}
	
	$badge_info = Get_badge_info(\@var_info);
	print $out "$badge_info	", join "\t", @var_info;
	
	print  $out "\n";
	
	#Call_ClinSig($badge_info);

}

close $out;

#################################################################
################  Sub-routines ################
	####### Get Badge Info  #########
#################################################################


sub Get_badge_info () {
	my ($arrRef) = @_; #pull in the array (@var_info)
    my @arr = @{$arrRef}; # conver array
    my @sig; 
    my %iSig;  ## counting the num of cases
    my @noBadge;
    my %inoBadge;
    my %ReviewPan;
    my $Var_call = "NA	NA";
    my $Var_call_exp = "NA	NA"; #expert call - var set up
    my $Alt_call = "NA	NA"; # alternant P-LP B-LB call (based on secondary rules (Rule2 -see later))
    my $P_LP_evid = "NA";  ## var call - set up
   	my %iSig1;
    my $other_badge_flag = "NA";
    my $other_exp_flag = "NA";
    my $other_Nonbadge_flag = "NA";
	my $return = "BadgeInfo-";
	my $return1  = "NonBadgeInfo-";
	my $return2  = "ReviewPanel-";
	my $rh = "NA"; ## flag for other affects
	my $rh1 = "NA"; ## flag for other affects
	my $rh2 = "NA"; ## flag for other affects
	my @PLP_dis_array;
	my @BLB_dis_array;
	my @vus_dis_array;
	my $res_conflict = "NA	NA";
	
	foreach $ln (@arr) {	
		@sig = split '\|', $ln;
		@di = split ':', $sig[4];   ### disease lookup
		next if $sig[0] < 20099999 and $sig[0] > 1;
		next if $sig[0] =~ /provided|Provided|not provided/;
		### Collect badge lab info
		#if ($ln =~ /Ambry|ARUP|Athena Diagnostics|Children's Hospital of Chicago|Robert H. Lurie|Children's Mercy Hospital|Color|GeneDx|GeneKor|Illumina|Integrated Genetics\/Laboratory Corporation of America|Invitae|Myriad Women's Health|Partners Laboratory for Molecular Medicine|Quest Diagnostics|University of Chicago/
		if ($ln =~ /$badgeList/
		) {	

					
			my $set;
			$iSig{$sig[1]}{$cnt}++;
			$iSig1{$sig[1]}{$di[0]}{$sig[0]} = $sig[0]; # get dates and disease
			if ($ln =~ /Pathogenic|Likely pathogenic|Likely_pathogenic/) {				
				$set = "$di[0]	$sig[0]	$sig[1]";  ## take up "Disease	Date	significance"				
				push @PLP_dis_array, $set; #put set inf into an array
			}
			elsif ($ln =~ /Benign|Likely benign|Likely_benign/) {				
				$set = "$di[0]	$sig[0]	$sig[1]";  ## take up "Disease	Date	significance"
				push @BLB_dis_array, $set;#put set inf into an array
			}
			
			elsif ($ln =~ /Uncertain significance/) {				
				$set = "$di[0]	$sig[0]	$sig[1]";  ## take up "Disease	Date	significance"
				push @vus_dis_array, $set;#put set inf into an array
			}

		}
		elsif ($ln =~ /expert panel|practice guideline/ and $_ !~ /ClinGen-approved/) { ### removing ClinGen-approved panels except for the Cardiomyopathy
				$ReviewPan{$sig[1]}{$cnt}++; 
			}
		
		else {  ### Collect non-badge lab info
		 	$inoBadge{$sig[1]}{$cnt}++;
		 }	
		 
		 ### Collect expert panel review information
	
	}
	
	
	#### Parse the badge lab information
	foreach $z ( %iSig) {
		if ($iSig{$z}{$cnt} > 0) {
			$return = $return."$z=$iSig{$z}{$cnt}:";		
			}
		}
	### Parse non-badge lab info
	foreach $y ( %inoBadge) {
		if ($inoBadge{$y}{$cnt} > 0) {
			$return1 = $return1."$y=$inoBadge{$y}{$cnt}:";	
			}	
		}
		
		### Parse expert panel review information
	foreach $x ( %ReviewPan) {
		if ($ReviewPan{$x}{$cnt} > 0) {
			$return2 = $return2."$x=$ReviewPan{$x}{$cnt}:";	
			}	
		}	

#################################################################
######################## Set flags for affects by expert panel
#################################################################
			
	foreach $rh (keys %ReviewPan) {
		if ($rh =~ /drug response|risk factor|protective|Affects|association/) {
			$other_exp_flag = "$rh";
			#print "$k	Found ReviewPanel 	$rh	\n";
		}
	}	
	
#################################################################	
######################## Set flags for affects by Badge_lab panel
#################################################################
		
	foreach $rh1 (keys %iSig1) {
		if ($rh1 =~ /drug response|risk factor|protective|Affects|association/) {
			$other_badge_flag = "$rh1";
			#print "$k	Found BadgeLab 		$rh1	\n";
		}
	}	


#################################################################	
######################## Set flags for affects by Badge_lab panel
#################################################################
		
	foreach $rh2 (keys %inoBadge) {
		if ($rh2 =~ /drug response|risk factor|protective|Affects|association/) {
			$other_Nonbadge_flag = "$rh2";
			#print "$k	Found NonBadgeLab 		$rh2	\n";
		}
	}	

	
	
#################################################################	
######################## count badge info
#################################################################
	
	my ($pathog, $pathog1, $pathog2, $benign, $VUS, $other, $benign1, $benign2, $VUS1, $pDiff, $pVUSdiff) = (0) x 11;
	$pathog = 	$iSig{Pathogenic}{$cnt} + $iSig{"Likely pathogenic"}{$cnt} + $iSig{"Pathogenic/Likely_pathogenic"}{$cnt};
	$pathog1 = $ReviewPan{Pathogenic}{$cnt} + $ReviewPan{"Likely pathogenic"}{$cnt} + $ReviewPan{"Pathogenic/Likely_pathogenic"}{$cnt};
	$pathog2 = $inoBadge{Pathogenic}{$cnt} + $inoBadge{"Likely pathogenic"}{$cnt} + $inoBadge{"Pathogenic/Likely pathogenic"}{$cnt} + $inoBadge{"Likely_pathogenic"}{$cnt};
	
	$benign = $iSig{Benign}{$cnt} + $iSig{"Likely_benign"}{$cnt} + $iSig{"Likely benign"}{$cnt} + $iSig{"Benign/Likely benign"}{$cnt} ;
	$benign1 = $ReviewPan{Benign}{$cnt} + $ReviewPan{"Likely_benign"}{$cnt} + $ReviewPan{"Benign/Likely benign"}{$cnt} ;
	$benign2 = $inoBadge{Benign}{$cnt} + $inoBadge{"Likely_benign"}{$cnt} + $inoBadge{"Benign/Likely benign"}{$cnt} ;

	$VUS = $iSig{"Uncertain significance"}{$cnt};
	$VUS1 = $inoBadge{"Uncertain significance"}{$cnt};
	$pDiff = $pathog2 - $benign2;
	$pVUSdiff = $pathog2 - $VUS1;
	
	my (@P_disease, @LP_disease);
	my %BLB;
	my $test = 0;
	my $Comp_P = 0;
	my $Comp_LP = 0;
	my $calls;


#################################################################	
	## Calling P or LP evidence. See later. 
#################################################################
		
	if ($pathog > 0 or $pathog1 > 0 or $pathog2 > 0 ) { $P_LP_evid = "Evidence_of_P/LP"; }
	
	else { $P_LP_evid = "No_evidence_of_P/LP"; }
	
	
#################################################################
##########   Call pathogenic / likely pathogenic	- Review Panel
#################################################################
	
	
	if ( $ReviewPan{"Pathogenic"}{$cnt} > 0 and  $ReviewPan{Pathogenic}{$cnt} => $ReviewPan{"Likely pathogenic"}{$cnt} ){
		$Var_call_exp = "Pathogenic	ReviewPanel_lab_P/LP";
		#print "Pathogenic	ReviewPanel_lab_P/LP	";
	}
	elsif ( $ReviewPan{"Pathogenic"}{$cnt} > 0  ){
		$Var_call_exp = "Pathogenic	ReviewPanel_lab_P/LP";
		#print "Pathogenic	ReviewPanel_lab_P/LP	";
	}
	
	elsif ( $ReviewPan{"Uncertain significance"}{$cnt} > 0 ){
		$Var_call_exp = "Uncertain significance	ReviewPanel_lab_VUS";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}
	
	elsif ( $ReviewPan{"Likely pathogenic"}{$cnt} > 0 ){
		$Var_call_exp = "Likely pathogenic	ReviewPanel_lab_P/LP";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}
	elsif ( $ReviewPan{"Likely_pathogenic"}{$cnt} > 0 ){
		$Var_call_exp = "Likely pathogenic	ReviewPanel_lab_P/LP";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}

	
	elsif (  $ReviewPan{"Benign"}{$cnt} > 0  ){
		$Var_call_exp = "Benign	ReviewPanelFinding_Benign";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}
	
	elsif ( $ReviewPan{"Likely benign"}{$cnt} > 0   ){
		$Var_call_exp = "Likely benign	ReviewPanelFinding_LB";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}
	
		elsif ( $ReviewPan{"Likely_benign"}{$cnt} > 0   ){
		$Var_call_exp = "Likely benign	ReviewPanelFinding_LB";
		#print "Likely_pathogenic	ReviewPanel_lab_P/LP	";
	}
	
	

#################################################################
##########   Call pathogenic / likely pathogenic	consensus (majority) - Badge labs
#################################################################

	elsif ($pathog > $benign and  $pathog > $VUS ) {  ### call P or LP - greater then VUS and B-LB
		
		if ($iSig{Pathogenic}{$cnt} > 0) {
			$Var_call = "Pathogenic	Badge-lab consenses P/LP";
		}
		else {$Var_call = "Likely pathogenic	Badge-lab consensus is P/LP"; }
		
	} # P-LP greater then VUS or BLB, and both VUS and BLB == 0
		
	elsif ( $pathog < $benign and $VUS < $benign  ) { ### call B or LB
		#print "Found Found	$k	No Conflict but more B/LB\n";
		if ($iSig{Benign}{$cnt} > 0) {
			$Var_call = "Benign	Badge-lab consensus is B/LB";
		}
		else {$Var_call = "Likely benign	Badge-lab consensus is B/LB"; }

	}  
	

#################################################################
##########   New rule for VUS, if P-LP is 0 use first (VUS consensus), if P-LP greater then 0 use rule 2
##########   This rule is separated from if VUS and P-LP are equal
##########   For rule2 use $Alt_call variable rather then $Var_call
##########    As per Sharon, Jung and Doug's request 12-3-2020 - We keep the original call for rule 2 but flag and provide an ALT call. 
#################################################################

	
	elsif ( $VUS > $pathog and $VUS > $benign and $pathog > 0.1 ) {## VUS rule 2
		my @rule2; #setting rule 2
		$Var_call = "Uncertain significance	Badge_lab consensus Uncertain significance - *Check secondary call";
		$Var_call1 = all_equal_conflict(\@PLP_dis_array,\@BLB_dis_array,\@vus_dis_array);
		@rule2 = split "\t", $Var_call1;
		$rule2[1] = "*Call based on having P-LP calls with VUS calls $rule2[1] - check before use - P-LP=$pathog B-LB=$benign VUS=$VUS";
		$Alt_call = join "\t", @rule2;
		
		#notify when rule is used
		#print "Invoking rule 2 $k	based on P-LP found Badge	vus=$VUS P-LP=$pathog	using all_equal_conflict\n";
		#print ">>>>>>>>>> $k	Calling $Var_call\n";

	} #### VUS > PLP and VUS > BLB
	
	elsif ( $VUS > $pathog and $VUS > $benign  ) { ## VUS rule 1 
	
		$Var_call = "Uncertain significance	Badge-lab consensus is Uncertain_significance";	
		
	} #### VUS > PLP and VUS > BLB


#################################################################
##########   Resolve conflicting 	- Badge labs B/LB and P/LP equal
#################################################################
	
	elsif ($pathog == $benign and $pathog > $VUS and $pathog > 0) { ### call PLP vs BLB conflict	
		#print "Check_disease	Check\n";
		$Var_call = Call_conflict(\@PLP_dis_array,\@BLB_dis_array);
		
	}
	
#################################################################
##########   Resolve conflicting 	- Badge labs VUS and P/LP equal
#################################################################

	
	elsif ($pathog == $VUS and $pathog > $benign and $pathog > 0) { ### call PLP vs VUS conflict

		$Var_call = Call_conflict(\@PLP_dis_array,\@vus_dis_array);
		#print "$k	Conflicts	$Var_call\n";		
	}
	

#################################################################
##########   Resolve conflicting 	- Badge labs VUS and B/LB equal
#################################################################
	
	elsif ($benign == $VUS  and $pathog < $benign and $benign > 0) { ### call BLB vs VUS conflict

		#print "$k	Check_disease	Check\n";

		$Var_call = Call_conflict(\@BLB_dis_array,\@vus_dis_array);

		#print "$k	Conflicts	$Var_call\n";
		
	}
	

#################################################################
##########   Resolve conflicting 	- Badge labs All equal - spacial case -- 


		elsif ($benign == $VUS  and $pathog == $benign and $benign > 0) { ### call BLB vs VUS conflict
			#$Var_call = "Conflicting_and_P	Check_disease";	

			$Var_call = all_equal_conflict(\@PLP_dis_array,\@BLB_dis_array,\@vus_dis_array);
		
	}


#################################################################
##########   ALT calls for nonBadge Lab	- Badge labs All equal - Use $Alt_call for these rather then $Var_call
#################################################################

	
	elsif ($pathog == 0 and $benign == 0 and $VUS == 0 and $pDiff > 1.9 and $pVUSdiff > 0.9  ){ # and $benign1 < 1.1 and $VUS1 < 1.1 and $pathog2 > 2.9
		if ($inoBadge{Pathogenic}{$cnt} > 1.1) {
			$Alt_call = "Pathogenic	*Non-badge labs consensus P/LP*. No badge-lab call";
		}
		else {$Alt_call = "Likely pathogenic	*Non-badge labs consensus P/LP*. No badge-lab call"; }
		
		#print "$k	Uncertain significance	Badge_labs_consensus_Uncertain_significance\n$Var_call_exp\n";
	}
	
	elsif ($pathog == 0 and $benign == 0 and $VUS == 0 and $pathog2 > 1.1 and $benign1 < 3 and $VUS1 < 1.1){
		if ($inoBadge{Benign}{$cnt} > 1.1) {
			$Alt_call = "Benign	*Non-badge labs consensus P/LP*. No badge-lab call";
		}
		else {$Alt_call = "Likely benign	*Non-badge labs consensus P/LP*. No badge-lab call"; }
		
		#print "$k	VUS	Badge_labs_all_were_Uncertain_significance\n$Var_call_exp\n";
	}
	
	elsif ($pathog == 0 and $benign == 0 and $VUS == 0 ){
		$Var_call = "No_info_available	There is no badge-lab nor Expert panel call available";
		#print "$k	VUS	Badge_labs_all_were_Uncertain_significance\n$Var_call_exp\n";
	}
	
	else {$Var_call = "No_info_available	There is no badge-lab, Expert panel or NonBadged-lab call available";  
		#print "$k	look_at_these\n";
		}
	
	

	if ($Var_call_exp ne "NA	NA" ) {
		$Var_call = $Var_call_exp;
		#print "$k	$Var_call_exp\n";	
	}
	if ($Alt_call eq "NA	NA" ) {
		$Alt_call = $Var_call;
		#print "$k	$Var_call_exp\n";	
	}

	if ($Var_call eq "NA	NA" ) {
		$Var_call = "No_info_available	There is no badge-lab, Expert panel or NonBadged-lab call available";
		#print "$k	$Var_call_exp\n";	
	}
	
	my $total_ret = "$P_LP_evid	$Var_call	$Alt_call	$other_exp_flag	$other_badge_flag	$other_Nonbadge_flag	$return	$return2	$return1";
	#print "$pathog,$benign,$VUS,$other	$Var_call	$P_LP_evid	$tot\n";
	#print "\n";
	return $total_ret;

}  # end of sub


#################################################################
########## Subroutine for Resolve conflicting 	- Badge labs - from line 342
#################################################################

sub Call_conflict  {
	my ($one_ref, $two_ref) = @_;
	my @P_conf = @{ $one_ref };       # dereferencing and copying each array
    my @B_conf = @{ $two_ref };
    
    my $return_Confl = "NA	NA";	
    my @l_conf;
    my %h_conf;
    my @date_P;
    my @date_B;
    
    
    $PPP = join '||', @P_conf; #outputting P-LP
    $BBB = join '||', @B_conf; #outputting B-LB
   
   
    foreach my $B (@B_conf) { # checking B-LB, getting dates
    	@l_conf = split "\t", $B;
    	$h_conf{$l_conf[0]}{$l_conf[1]} = $l_conf[1]; 
    	push @date_B, "$l_conf[1]	$l_conf[2]	$l_conf[0]";
    	
    	#print "Check Disease	$l_conf[0]	$l_conf[1]	$l_conf[2]\n";  	#checking on condition
    }
	
	foreach my $P (@P_conf) { #organizing/sorting the data - P-LP
    	@l_conf = split "\t", $P;
    	push @date_P, "$l_conf[1]	$l_conf[2]	$l_conf[0]";   		
	}

	##### Comparing diseases - if a disease for P-LP found that is not in others - call P-LP	
	foreach my $P (@P_conf) {
		@l_conf = split "\t", $P;
		
		if  (!exists $h_conf{$l_conf[0]} and $l_conf[0] !~ /CN517202|CN169374|not specified|Specified|not provided/) {
    		$return_Confl = "$l_conf[2]	Conflicting but disease $l_conf[0] was found in P-LP and not in B-LP or VUS";
    		#print "******** Call_conflict	$k	$l_conf[2]	Conflicting but disease $l_conf[0] was found in $l_conf[2]  P=$PPP	||	B=$BBB\n"; 
    		last;    		  		
   	 	}
   	 	
   	}

	if ($return_Confl =~ /NA/) { # if no resolution on disease
			#foreach $mm (@B_conf) {print "getInfoConf $mm\t";}
		#	my @sorted = sort {(split('	',$b))[1]<=>(split('	',$a))[1]} @array; #Example
		@date_Ps = sort {(split('	',$b))[1]<=>(split('	',$a))[1]} @P_conf;
		@date_Bs = sort {(split('	',$b))[1]<=>(split('	',$a))[1]} @B_conf;
		my ($dip, $dp,  $sp) = split "\t", $date_Ps[0];
		my ($dib,$db, $sb) = split "\t", $date_Bs[0];

		
		if ($dp > $db) {
			$return_Confl = "$sp	*Latest reported date call is $sp date $dp";
			#print "$k	$sp			||		$dp $sp $dip is later then $db $sb $dib	OtherEnd= $date_Ps[-1]\t";

			
		}
		elsif ($dp < $db) {
			$return_Confl = "$sb	*Latest reported date call is $sb date $db";
			#print "$k	$sb				||		$db $sb $dib is later then $dp $sp $dip	OtherEnd= $date_Ps[-1]\n";

		}
		
		else {
			#print "$k	Call_conflict	Can't resolve conflict	No clear resolution	Call_conflict\n"; 
			$return_Confl = "Can't resolve conflict	No clear resolution";
		}
		
	
	}

	return $return_Confl;

}  ### end of sub

#################################################################
##########   Resolve conflicting 	- Badge labs All equal - spacial case -line 380
#################################################################

sub all_equal_conflict {
		
	my ($one_ref, $two_ref, $three_ref) = @_;
	# dereferencing and copying each array
	## each has the following info 1.Disease	2.Date	3.significance
	my @P_conf = @{ $one_ref };    # P LP array    
    my @Ba_conf = @{ $two_ref };	# B LB array 
    my @V_conf = @{ $three_ref }; 	# VUS array 
    
    ### combining VUS and B-LB - since they are compared with P-LP
    ### Used only for disease comparison - if there is a P-LP disease not in others. 
    my @B_conf = (@Ba_conf, @V_conf); 
     
 
	my $return_Confl = "NA	NA";	# set the result to NA
    my @l_conf;
    my %h_conf;
    my @date_P;
    my @date_B;
    
    #print "$k	P-LP, VUS, B-LB are equal. Checking for resolution!!!!!!!!!!!!!!!!!!!!!\t"; #check print
    
    $PPP = join '||', @P_conf;
    $BBB = join '||', @B_conf;
   
    # $set = "$di[0]	$sig[0]	$sig[1]";  ## take up "Disease	Date	significance"

    foreach my $B (@B_conf) { #organizing the data - B-LB
    	@l_conf = split "\t", $B;
    	$h_conf{$l_conf[0]}{$l_conf[1]} = $l_conf[1]; 
    	push @date_B, "$l_conf[1]	$l_conf[2]	$l_conf[0]";
    	
    	#print "Check Disease	$l_conf[0]=disease	$l_conf[1]=date	$l_conf[2]=sig\n";  	
    }
	
	foreach my $P (@P_conf) { #organizing/sorting the data - P-LP
    	@l_conf = split "\t", $P;
    	push @date_P, "$l_conf[1]	$l_conf[2]	$l_conf[0]";   		
	}

	##### Comparing diseases - if a disease for P-LP found that is not in others - call P-LP	
	foreach my $P (@P_conf) {
		@l_conf = split "\t", $P;
		
		if  (!exists $h_conf{$l_conf[0]} and $l_conf[0] !~ /CN517202|CN169374|not specified|Specified|not provided/) {
    		$return_Confl = "$l_conf[2]	Conflicting but disease $l_conf[0] was found with $l_conf[2]";
    		#print "********$k	$l_conf[2]	Conflicting but disease $l_conf[0] was found with $l_conf[2] or VUS  P=$PPP	||	B=$BBB\n"; 
    		last;    		  		
   	 	}
   	 	
   	}
    	   
	### if still not resolved - go by date, get the latest date for P-LP and B-LB - if still the same goto cant reolve
	if ($return_Confl =~ /NA/) {
		#foreach $mm (@B_conf) {print "getInfoConf $mm\t";}
		
		@date_Ps = sort {(split('	',$b))[1]<=>(split('	',$a))[1]} @P_conf;
		@date_Bs = sort {(split('	',$b))[1]<=>(split('	',$a))[1]} @B_conf;
		my ($dip, $dp, $sp) = split "\t", $date_Ps[0];
		my ($dib, $db, $sb) = split "\t", $date_Bs[0];
		
		if ($dp > $db) {
			$return_Confl = "$sp	*Latest reported date is $dp for call $sp";
			#print "$k	$sp			||		$dp $sp $dip is later then $db $sb $dib	OtherEnd= $date_Ps[-1]\n";
			
		}
		elsif ($dp < $db) {
			$return_Confl = "$sb	*Latest reported date is $db for call $sb";
			#print "$k	$sb				||		$db $sb $dib is later then $dp $sp $dip OtherEnd= $date_Ps[-1]\n";
			
		}
		
		else {
			#print "$k	Can't resolve conflict	No clear resolution	all_equal_conflict\n"; 
			$return_Confl = "Can't resolve conflict	No clear resolution";
		}
		
	
	}

	return $return_Confl;
}











