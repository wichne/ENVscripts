#!/usr/bin/perl

print LOG "Hello. I am the GenomeDaemon. I have found a new fasta file to process called ${FASTA_FILE}.
    First, let's look at some information about this file.";

$INPUT_FASTA_STATS = &fasta_stats($FASTA_FILE);

print LOG "The file contains $INPUT_FASTA_STATS->{num_records} sequence records.\n";
print LOG "The total sequence length is $INPUT_FASTA_STATS->{tot_length} bp.\n";
print LOG "The shortest sequence is $INPUT_FASTA_STATS->{shortest} bp.\n";
print LOG "The shortest sequence is $INPUT_FASTA_STATS->{shortest} bp.\n";
print LOG "The longest sequence is $INPUT_FASTA_STATS->{longest} bp.\n";
print LOG "Records >1kb comprise $INPUT_FASTA_STATS->{short_percent} % of the total length.\n";
print LOG "The characters represented in the sequence are ['" . join("', '", @{$INPUT_FASTA_STATS->{alphabet}}) . "'].\n";
print LOG "This looks like $INPUT_FASTA_STATS->{dbtype} sequence.\n";

print LOG "The %G+C for the total sequence is $INPUT_FASTA_STATS->{GC} %.\n";

# bunch of ifs...

print LOG "This looks like a genome sequence. Let's begin.\n";
print LOG "First we will run PROKKA, since it nicely combines many relevant gene-finding steps and a baseline annotation\n";
# Run PROKKA and load
print LOG "PROKKA identified:\ntype\tcount\tsuccessfully loaded in to $DB\n";
foreach $feature_type (@PROKKA_feature_types) {
    print LOG "$feature_type\t$PROKKA_OUTPUT->{$feature_type}\t$PROKKA_LOAD_success->{$feature_type}\n";
}
# Evidence
