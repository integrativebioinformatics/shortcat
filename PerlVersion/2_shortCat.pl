#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# NOTE: To Strand specific analysis, to modify command found in line 41

my $usage = "
Usage: shortCat OPTION:

 -i, --input	  sam file to calculate reads count matrix
 -c, --conditions comma separated list of conditions
 -o, --output	  output file name [Default: same input with extension .count]
 -s, --strand     strand specific data (yes, no) [Default: no]

Example:
	shortCat -i reads_mapped.sam -c KO1,KO2,WT1,WT2
*where KO1,KO2,WT1,WT2 are the read prefixes in fastq headers (ex: WT1_Numbers)

";

my ($help, $filein, $fileout, $conditions );
my $strand_specific = "no";
if ( @ARGV < 1 or !GetOptions('help|h' => \$help, 'input|i=s' => \$filein, 'output|o=s' => \$fileout, 'conditions|c=s' => \$conditions, 'strand|s=s' => \$strand_specific)
          or defined $help ) {
        print "Unknown option: @_\n" if ( @_ );
 	print $usage;
	exit;
}
if(!defined($filein) || !defined($conditions)) {
        print "Unknown option: @_\n" if ( @_ );
 	print $usage;
	exit;
}
if ($strand_specific !~ /yes|no/){
	print "Unknown option: Strand specification\n";
	print $usage;
	exit;
}
open(ACE,$filein) or die $!;
close(ACE);
my @filename = split(/\./,$filein);
if(!defined($fileout)) {
	$fileout = $filename[0].".count";
}
open(PRI,">$fileout");
my @condition = split(",",$conditions);

system("sam2bed --max-mem 70G < $filein >$filename[0].bed");

if ($strand_specific eq "yes"){
	system("mergeBed -s -i $filename[0].bed -c 4,6 -o distinct -delim \";\" >$filename[0]_merge.bed");
} else {
	system("mergeBed -i $filename[0].bed -c 4,6 -o distinct -delim \";\" >$filename[0]_merge.bed");
}

#system("mergeBed -i $filename[0].bed -nms >$filename[0]_merge.bed");

open(ACE,"$filename[0]_merge.bed") or die $!;

my $count = 1;
#PreWt,PreTg,SinWt,SinTg
print PRI "Chrom\tStart\tEnd\tContig\tLength\tstrand\t#reads";
foreach my $condi (@condition) {
	print PRI "\t$condi";
}
print PRI "\n";
my $indicador = 0;
while(<ACE>){
    chomp;
    my @col = split(/\t/,$_);
    my @infoAf = split(/\;/,$col[3]);
    my %counts = ();
    foreach my $condi (@condition) {
	$counts{$condi} = 0;		#cambiar a 1 para sumarle 1 a todas las cuentas por defecto
    }
    my $cantidad = @infoAf;
#   if ($cantidad > 1) {
	foreach my $ests (@infoAf) {
		$indicador = 0;
		foreach my $condi (@condition) {
			if($ests =~ /^$condi/) {
				$counts{$condi}++;
				$indicador = 1;
			}
		}
		if ($indicador == 0) {
			die "Tipo no encontrado : $_\tline: $count";
		}
	} 
	my $largo = $col[2] - $col[1];
	print PRI "$col[0]\t$col[1]\t$col[2]\tSmallRNA_$count\t$largo\t$col[4]\t$cantidad";
	foreach my $condi (@condition) {	
		print PRI "\t$counts{$condi}";
	}
	print PRI "\n";

#   }
    $count++;
}
close(PRI);
