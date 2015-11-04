#!/usr/bin/perl

use DBI;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib "/home/wcnelson/devel";
use ENV;

my $args = {};
&getopts('f:D:', $args);

my $filename = $args->{'f'} or die "Need to provide filename with -f\n";

my $dbh = DBI->connect("dbi:mysql:server=localhost;database=$args->{D}",
		       $ENV{USER}, "RenMan", {'AutoCommit'=> 0});

#open file
open my $in, $filename || die "Can't open $filename\n";

my %data;
my $seq_count;
while (my $line = <$in>) {
    chomp $line;
    my ($key, $value) = split/\t/, $line;
    if ($key eq "ca_contig_id") { $data{seq_accession} = $value }
    elsif ($key eq "com_name") { $data{description} = $value }
    elsif ($key eq "seq#") { $data{seq_count} = $value }
    elsif ($key eq "seq_name") {
	my %subseq;
	$seq_count++;
	$subseq{$key} = $value;
	until ($line eq "" || $line eq "|") {
	    chomp($line = <$in>);
	    my ($key, $value) = split/\t/, $line;
	    $subseq{$key} = $value;
	}

	# this is because we only loaded the clear range, so the tasm coords
	# are not right.
#	my ($sub_min, $sub_max) = sort {$a<=>$b} ($subseq{seq_lend}, $subseq{seq_rend});
	my $revcomp = $subseq{seq_lend} > $subseq{seq_rend} ? 1 : 0;
	my $query = "SELECT a1.sequence_id, $subseq{asm_lend}, $subseq{asm_rend},"
	    . " a2.sequence_id, 1, s.sequence_length, $revcomp"
	    . " FROM sequence_accessions a1, sequence_accessions a2, sequences s"
	    . " WHERE a1.seq_accession='$data{seq_accession}'"
	    . " AND a2.seq_accession='$subseq{seq_name}'"
	    . " AND s.sequence_id=a2.sequence_id;\n"
	    ;
	print $query;
    }
    
    if ($key eq "|") {
	%data = ();
	$seq_count = 0;
    }
}
