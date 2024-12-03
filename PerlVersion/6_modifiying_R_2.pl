## THIS SCRIPT MODIFY THE R SCRIPT OF DESEQ2 BASED ON ANALISYS OF SHORTCAT

## USAGE: perl $0 deseq2_script RUTE_COUNTS RUTE_CONDITION CONDITION_FILE

open ($fh_1,">","No_1.R");
open ($fh_2,">","No_2.R");

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

## --------------------------------------------------------------------------------------------------------------

#library("DESeq2")
#cts <- as.matrix(read.csv("all_filter_cutoff_5_toDESEQ2.tab",sep="\t",row.names="gene_id"))
#coldata <- read.csv("CONDITION.csv",sep="\t", row.names=1)
#coldata <- coldata[,c("condition","type")]
#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = coldata,
#                              design = ~ condition)
#dds <- DESeq(dds)
#table_counts_normalized <- counts(dds, normalized=TRUE)
#write.csv(as.data.frame(table_counts_normalized), file="DESEQ2_normalized_sRNAs.csv")

## --------------------------------------------------------------------------------------------------------------

#library("DESeq2")
#cts <- as.matrix(read.csv("all_filter_cutoff_5_toDESEQ2.tab",sep="\t",row.names="gene_id"))
#coldata <- read.csv("CONDITION.csv",sep="\t", row.names=1)
#coldata <- coldata[,c("condition","type")]
#bases <- as.matrix(read.table("sRNA_length.tab",sep="\t",row.names=1, header = T))
#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = coldata,
#                              design = ~ condition)
#mcols(dds)$basepairs <- bases[, "Length"]
#fpkm <- fpkm(dds, robust = TRUE)
#write.csv(as.data.frame(fpkm), file="fpkm_sRNAs.csv")

#library("DESeq2")
#cts <- as.matrix(read.csv("COUNTS",sep="\t",row.names="gene_id"))
#coldata <- read.csv("CONDITIONS",sep="\t", row.names=1)
#coldata <- coldata[,c("condition","type")]
#bases <- as.matrix(read.table("z_sRNAs_length.tab",sep="\t",row.names=1, header = T))
#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = coldata,
#                              design = ~ condition)
#mcols(dds)$basepairs <- bases[, "Length"]
#fpkm <- fpkm(dds, robust = TRUE)
#write.csv(as.data.frame(fpkm), file="FPKM_normalized_sRNAs.csv")
