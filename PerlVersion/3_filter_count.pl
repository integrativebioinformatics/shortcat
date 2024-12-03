## USAGE: perl 3_filter_count.pl TAB_FILE COUNT_FILE CUTOFF_READS

##	TAB_FILE
##		File of two columns having as column 1 whatever (ex: name file) and column 2 the ID used for SHORTCAT

##	COUNT_FILE
##		File produced by SHORTCAT

##	CUTOFF_READS
##		Number of reads necessaries to consider the sRNA in a sample

##EXAMPLE USAGE:
##	perl 3_filter_count.pl ../z_script1_ruteCATfiles_R1.tab all_prueba.count 5



## This script to generate the filter of all.count file produced by SHORTCAT
## 	Requires that ID of samples for SHORTCAT not contain NUMBERS, with exception of numeration of replicates, because to select the samples that are replicates, the command will remove numbers from ID samples for SHORCAT. 

##	EXAMPLES OF VALIDS ID (column 2):

##		Name_sample_1	CTRL1		Script process end with CTRL
## 		Name_sample_2	CTRL2		Script process end with CTRL

##		Name_sample_1	CTRL_1		Script process end with CTRL_
## 		Name_sample_2	CTRL_2		Script process end with CTRL_

##		Name_sample_1	CTRL12		Script process end with CTRL
## 		Name_sample_2	CTRL13		Script process end with CTRL


##	EXAMPLES OF NO VALIDS ID (column 2):

##		Name_sample_1	C1TRL1		Script process end with CTRL (this will not found in the matrix counts)
## 		Name_sample_2	C1TRL2		Script process end with CTRL (this will not found in the matrix counts)

##		Name_sample_1	5CTRL_1		Script process end with CTRL_ (this will not found in the matrix counts)
## 		Name_sample_2	5CTRL_2		Script process end with CTRL_ (this will not found in the matrix counts)

##		Name_sample_1	57cr82_1	Script process end with cr_ (this will not found in the matrix counts)
## 		Name_sample_2	57cr82_2	Script process end with cr_ (this will not found in the matrix counts)

## NOTE:
##	The script report only the formed contigs that in
##	AT LEAST ONE CONDITION (i.e. all replicates of one condition)
##	THE 50% OF SAMPLES HAVE MORE THAN "N" READS FOR THAT CONTIG. ("N" is defined in ARGV[2])




## ----------------------------------------------- PART I ----------------------------------------------- ##
## 		Saving the column number of each sample respect to ID of replicates 




## -- RECOVER ID REPLICATES-- ##
my @id = `cut -f 2 $ARGV[0] | sed -e 's/[0-9]\\+//g' | sort -u`;
chomp @id;

my $head_counts = `head -1 $ARGV[1]`;
chomp $head_counts;

print "$head_counts\n";

my @div_head = split ("\t", $head_counts);


## -- MAKE HASH OF HASH TO SAVE NUMBER OF COLUMN OF EACH REPLICATE-- ##
my %ID_1 = ();

foreach my $each_id (@id){

	my $index_1 = 0;

	foreach my $each_col_head (@div_head){
		if ($each_col_head =~ $each_id){
			$ID_1{$each_id}{$each_col_head} = $index_1; ## Here save the column number ex: "$ID_1{CRL}{CRL1} = 7"
		}
		$index_1++;
	}
}


## ----------------------------------------------- PART II ----------------------------------------------- ##
## Generating the parse based on a specific number of reads on each sample




## -- OPEN COUNT FILE -- ##
open (IN,"<",$ARGV[1]);

<IN>;

while (<IN>){
	chomp $_;

	my $count_final = 0; # Count the number of condition with >= 50% of samples with at least "N" reads in their contig formed

	my @cols_line = split ("\t", $_);

	foreach my $each_id_first (@id){
		my $number_samples_id = 0; # HERE the number of samples
		my $count_of_replicates = 0; # HERE the number of replicates with more than "N" reads
		foreach my $values (values %{$ID_1{$each_id_first}}){
			$number_samples_id++;
			if ($cols_line[$values] >= $ARGV[2]){
				$count_of_replicates++;
			}
		}

		my $percent = ($count_of_replicates * 100) / $number_samples_id; # HERE percent of samples

		if ($percent >= 50 ){
			$count_final++;
		}

	}

	if ($count_final >= 1){
		my $recover_line = join ("\t", @cols_line);
		print "$recover_line\n";
	}


}

close (IN);

































