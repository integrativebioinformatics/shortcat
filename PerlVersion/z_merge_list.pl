#!/usr/bin/perl
use strict;
use warnings;
if(@ARGV !=2) {
	print "\nUsage: merge_list complete_file info_file\n\n";
	exit;
}
my $DE = $ARGV[0];
my $info_file = $ARGV[1];

open(ACE,$info_file) or die $!;
my %information = ();
while(<ACE>){
	chomp;
	my @fields = split("\t",$_);
	my $contig = shift(@fields);
	my $info = join("\t",@fields);
	$information{$contig} = $info;		
}
close(ACE);
open(ACE2,$DE);
while(<ACE2>){
	chomp;
	my @fields = split("\t",$_);
	if($information{$fields[0]}) {
		print "$_\t$information{$fields[0]}\n";
	} else {
		print "$_\tNO_MERGE\n";
	}
}
close(ACE2);
