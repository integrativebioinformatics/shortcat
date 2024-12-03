open (IN,"<",$ARGV[0]);

while (<IN>){
	chomp $_;

	my ($sRNA_id, $annot) = split ("\t", $_);

	if ($annot =~ /\|/){
		my %remove_repeats;

		my @div_annot = split (/\|/, $annot);

		foreach my $each_annot (@div_annot){
			$remove_repeats{$each_annot} = 0;
		}

		my @recover_no_repeats = keys (%remove_repeats);

		my $join_all = join ("|", @recover_no_repeats);

		print "$sRNA_id\t$join_all\n";
	} else {
		print "$sRNA_id\t$annot\n";
	}
}
