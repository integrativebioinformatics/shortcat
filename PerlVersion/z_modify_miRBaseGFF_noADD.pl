if (@ARGV < 1){
	die "USAGE: perl $0 GFF_miRBase\n";
}

my @lines_gff_miRBase = `cat $ARGV[0]`;
chomp @lines_gff_miRBase;

foreach my $each_line_gff (@lines_gff_miRBase){

	if ($each_line_gff !~ /^#/){
		my @cols_line = split ("\t", $each_line_gff);
		my $annot = $cols_line[-1];

		my @div_annot = split ("\;", $annot);

		$cols_line[-1] = "NO_NAME";
			
		foreach my $each_annot (@div_annot){
			if ($each_annot =~ /Name=/i){
				$each_annot =~ s/Name=//g;
				$cols_line[-1] = $each_annot;
			}
		}

		my $recover_line = join ("\t", @cols_line);
		print "$recover_line\n";
	} else {
		print "$each_line_gff\n";
	}


}
