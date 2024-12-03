# USED AS:
# cat *fastqc.html | perl 2_adapter_seq.pl - | cut -f 1 -d "(" | sort -u | grep -v "Hit" |awk -F "\t" '{gsub(/ /,"_",$2); gsub(/,/,"",$2);print ">" $2 "\n" $1}' > 1_adapters_chipseq.fna

@line_adapter = `grep "Overrepresented sequences</h2><table><thead><tr><th>Sequence</th><th>Count</th><th>Percentage</th><th>Possible Source</th></tr></thead><tbody><tr><td>" $ARGV[0]`;
chomp @line_adapter;

foreach $each_line (@line_adapter){
	@div_line = split ("Overrepresented sequences</h2><table><thead><tr><th>Sequence</th><th>Count</th><th>Percentage</th><th>Possible Source</th></tr></thead><tbody><tr><td>", $each_line);
	$first_recover = pop (@div_line);
	@div_line_2 = split (/[bp\)|Hit]<\/td><\/tr><tr><td>/, $first_recover);
	$number_strings = scalar (@div_line_2);
	if ($number_strings == 1){
		foreach $each_element (@div_line_2){
			@div_element = split (/<\/td><td>/, $each_element);
			$sequence = $div_element[0];

			@div_rest_element = split (/<\/td><\/tr><\/tbody>/, $div_element[3]);
			
			$identification = $div_rest_element[0];
			$identification =~ s/No Hi/No Hit/g;
#			$identification =~ s/\(|\)//g;
			print "$sequence\t$identification\n";
		}
	} else {
		foreach $each_element (@div_line_2){
			@div_element = split (/<\/td><td>/, $each_element);
			$sequence = $div_element[0];
			
			@div_rest_element = split (/<\/td><\/tr><\/tbody>/, $div_element[3]);
			
			$identification = $div_rest_element[0];
			$identification =~ s/No Hi/No Hit/g;
#			$identification =~ s/\(|\)//g;
			print "$sequence\t$identification\n";
		}
	}
}
