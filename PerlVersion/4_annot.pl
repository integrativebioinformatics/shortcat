use Getopt::Long qw(GetOptions);
#use diagnostics;

if (@ARGV < 2){
die "\nUSAGE: perl $0 [options] FILE_ANNOT_1 FILE_ANNOT_2\n
		Options
		-f	Minimum overlap required as a fraction of A (default 50%). 
		-r	Require that the fraction overlap be reciprocal for A AND B.
		-s	Require same strand of both features\n\n";	
}




## ----- DEFINING THE OPTIONS

my $the_f_option = "0.50";
my $the_r_option;
my $the_s_option;

GetOptions(
	'from=f' => \$the_f_option,
	'r' => \$the_r_option,
	's' => \$the_s_option,
) or die "Not valid option\n";


## ----- Making Intersect

if ($the_r_option){
	if ($the_s_option){
		our @lines_intersect = `intersectBed -wao -s -r -f $the_f_option -a $ARGV[0] -b $ARGV[1]`;
		chomp @lines_intersect;

#		print "Command used was: intersectBed -wao -s -r -f $the_f_option -a $ARGV[0] -b $ARGV[1]\n"
	} else {
		our @lines_intersect = `intersectBed -wao -r -f $the_f_option -a $ARGV[0] -b $ARGV[1]`;
		chomp @lines_intersect;

#		print "Command used was: intersectBed -wao -r -f $the_f_option -a $ARGV[0] -b $ARGV[1]\n"
	}
} else {
	if ($the_s_option){
		our @lines_intersect = `intersectBed -wao -s -f $the_f_option -a $ARGV[0] -b $ARGV[1]`;
		chomp @lines_intersect;

#		print "intersectBed -wao -s -f $the_f_option -a $ARGV[0] -b $ARGV[1]\n"
	} else {
		our @lines_intersect = `intersectBed -wao -f $the_f_option -a $ARGV[0] -b $ARGV[1]`;
		chomp @lines_intersect;

#		print "intersectBed -wao -f $the_f_option -a $ARGV[0] -b $ARGV[1]\n"
	}
}

my %sRNA_annot_gff = ();

foreach my $each_line_inter (@lines_intersect){

	my @cols_inter = split ("\t", $each_line_inter);

	my $sRNA_id = $cols_inter[3];
	my $gff_id = $cols_inter[-2];

	if ($cols_inter[-1] != 0){

		if (!defined $sRNA_annot_gff{$sRNA_id}){
			$sRNA_annot_gff{$sRNA_id} = "$gff_id";
		} else {
			if ($gff_id ne $sRNA_annot_gff{$sRNA_id}){
				$sRNA_annot_gff{$sRNA_id} = "$sRNA_annot_gff{$sRNA_id}" . "|" . "$gff_id";
			}	
		}

	} else {
	
		$sRNA_annot_gff{$sRNA_id} = "-"

	}
}

foreach my $each_key_annot (keys %sRNA_annot_gff){

	print "$each_key_annot\t$sRNA_annot_gff{$each_key_annot}\n";

}


























