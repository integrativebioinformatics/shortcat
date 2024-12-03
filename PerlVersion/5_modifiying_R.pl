## THIS SCRIPT MODIFY THE R SCRIPT OF DESEQ2 BASED ON ANALISYS OF SHORTCAT

## USAGE: perl $0 deseq2_script RUTE_COUNTS RUTE_CONDITION CONDITION_FILE

open ($fh_x,">","05-DiffExpression/No_x.R");
open ($fh_1,">","05-DiffExpression/No_1.R");
open ($fh_2,">","05-DiffExpression/No_2.R");

## ------- DIFF EXP

if (defined $ARGV[0]){
	my @lines_deseq_script = `cat $ARGV[0]`;
	chomp @lines_deseq_script ;

	my $rute_file_counts = $ARGV[1];
	my $rute_file_condition = $ARGV[2];

	my @div_cts = split ("HERE_RUTE_COUNTS", $lines_deseq_script[1]);
	my $recover_cts = join ("$rute_file_counts", @div_cts);
	$lines_deseq_script[1] = $recover_cts;

	my @div_coldata = split ("HERE_RUTE_CONDITION", $lines_deseq_script[2]);
	my $recover_coldata = join ("$rute_file_condition", @div_coldata);
	$lines_deseq_script[2] = $recover_coldata;


	my @id = `cut -f 2 05-DiffExpression/$ARGV[2] | awk 'NR > 1' | sort -u | awk '{print "\\"" \$1 "\\""}'`;
	chomp @id;

	my @id_2 = `cut -f 2 05-DiffExpression/$ARGV[2] | awk 'NR > 1' | sort -u`;
	chomp @id_2;

	my $the_levels = join (",", @id);


	my @dds_condition = split ("HERE_THE_ID", $lines_deseq_script[7]);
	my $recover_dds_condition = join ("$the_levels", @dds_condition);
	$lines_deseq_script[7] = $recover_dds_condition;

	my $part_I_deseq2 = join ("\n", @lines_deseq_script);

	print $fh_x "$part_I_deseq2\n";




	my $size_array = scalar (@id_2);
	my $index_for_compare_outside = 1;

	foreach my $each_id_2 (@id_2){

		if ($index_for_compare_outside < $size_array){

			my $index_for_compare_inside = $index_for_compare_outside;

			for (my $i = 1; $i < $size_array; $i++){
				if (defined $id_2[$index_for_compare_inside]){
					print $fh_x "res <- results(dds, contrast=c(\"condition\",\"$each_id_2\",\"$id_2[$index_for_compare_inside]\"))\n";
					print $fh_x "write.csv(as.data.frame(res), file=\"$each_id_2\_vs_$id_2[$index_for_compare_inside].csv\")\n";
					$index_for_compare_inside++;
				}
			}

			$index_for_compare_outside++;

		}
	}
}


## ----------------------------------------------------------------------------------------------------

## ------- FPKM

if (defined $ARGV[0]){

	my @lines_deseq_script = `cat $ARGV[0]`;
	chomp @lines_deseq_script ;

	my $nothing_1 = pop (@lines_deseq_script);
	my $nothing_2 = pop (@lines_deseq_script);

	my $rute_file_counts = $ARGV[1];
	my $rute_file_condition = $ARGV[2];

	my @div_cts = split ("HERE_RUTE_COUNTS", $lines_deseq_script[1]);
	my $recover_cts = join ("$rute_file_counts", @div_cts);
	$lines_deseq_script[1] = $recover_cts;

	my @div_coldata = split ("HERE_RUTE_CONDITION", $lines_deseq_script[2]);
	my $recover_coldata = join ("$rute_file_condition", @div_coldata);
	$lines_deseq_script[2] = $recover_coldata;

	print $fh_1 "$lines_deseq_script[0]\n";
	print $fh_1 "$lines_deseq_script[1]\n";
	print $fh_1 "$lines_deseq_script[2]\n";
	print $fh_1 "$lines_deseq_script[3]\n";
	print $fh_1 "bases <- as.matrix(read.table(\"z_sRNAs_length.tab\",sep=\"\\t\",row.names=1, header = T))\n";
	print $fh_1 "$lines_deseq_script[4]\n";
	print $fh_1 "$lines_deseq_script[5]\n";
	print $fh_1 "$lines_deseq_script[6]\n";
	print $fh_1 "mcols(dds)\$basepairs <- bases[, \"Length\"]\n";
	print $fh_1 "fpkm <- fpkm(dds, robust = TRUE)\n";
	print $fh_1 "write.csv(as.data.frame(fpkm), file=\"FPKM_normalized_sRNAs.csv\")\n";

}

## ------- DESEQ2

if (defined $ARGV[0]){

	my @lines_deseq_script = `cat $ARGV[0]`;
	chomp @lines_deseq_script ;

	my $nothing_1 = pop (@lines_deseq_script);
	my $nothing_2 = pop (@lines_deseq_script);

	my $rute_file_counts = $ARGV[1];
	my $rute_file_condition = $ARGV[2];

	my @div_cts = split ("HERE_RUTE_COUNTS", $lines_deseq_script[1]);
	my $recover_cts = join ("$rute_file_counts", @div_cts);
	$lines_deseq_script[1] = $recover_cts;

	my @div_coldata = split ("HERE_RUTE_CONDITION", $lines_deseq_script[2]);
	my $recover_coldata = join ("$rute_file_condition", @div_coldata);
	$lines_deseq_script[2] = $recover_coldata;

	print $fh_2 "$lines_deseq_script[0]\n";
	print $fh_2 "$lines_deseq_script[1]\n";
	print $fh_2 "$lines_deseq_script[2]\n";
	print $fh_2 "$lines_deseq_script[3]\n";
	print $fh_2 "$lines_deseq_script[4]\n";
	print $fh_2 "$lines_deseq_script[5]\n";
	print $fh_2 "$lines_deseq_script[6]\n";
	print $fh_2 "dds <- DESeq(dds)\n";
	print $fh_2 "table_counts_normalized <- counts(dds, normalized=TRUE)\n";
	print $fh_2 "write.csv(as.data.frame(table_counts_normalized), file=\"DESEQ2_normalized_sRNAs.csv\")\n";

}

close ($fh_1);
close ($fh_2);










#library("DESeq2")
#cts <- as.matrix(read.csv("HERE_RUTE_COUNTS",sep="\t",row.names="gene_id"))
#coldata <- read.csv("HERE_RUTE_CONDITION",sep="\t", row.names=1)
#coldata <- coldata[,c("condition","type")]
#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = coldata,
#                              design = ~ condition)
#dds$condition <- factor(dds$condition, levels = c(HERE_THE_ID))
#dds <- DESeq(dds)



#res <- results(dds, contrast=c("condition","control","leve"))
#write.csv(as.data.frame(res), 
#          file="1_n1_m10_C_vs_L.csv")

#res <- results(dds, contrast=c("condition","control","grave"))
#write.csv(as.data.frame(res), file="2_n1_m10_C_vs_G.csv")

#res <- results(dds, contrast=c("condition","leve","grave"))
#write.csv(as.data.frame(res), 
#          file="3_n1_m10_L_vs_G.csv")
