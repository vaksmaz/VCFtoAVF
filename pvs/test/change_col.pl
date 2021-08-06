#!/usr/bin/perl




$filename = $ARGV[0];  # intervar

if ($filename =~ /gz$/){
	open($fh, "gunzip -c $filename |") or die "gunzip $filename: $!";
}
else {
	open $fh, $filename || die "Can't read file head.txt file '$filename' [$!]\n";
}


$header = <$fh>; chomp $header;

print "$header\n";

while (<$fh>) {
	chomp;
	@l = split "\t", $_;
	
	$l[2] = '.';
	
	print join "\t", @l, "\n";
	
}
