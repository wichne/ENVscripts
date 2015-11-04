#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

my $opt = {};
&getopts('D:i:n:o:', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'access');

my $OUT;
if ($outfile) {
    $OUT = Bio::SeqIO->new(-file => ">$outfile",
			   -format => 'genbank');
} else {
    $OUT = Bio::SeqIO->new(-fh => \*STDOUT,
			   -format => 'genbank');
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $seq_ids = &set_id_to_seq_ids($dbh, $setid);
foreach my $seqid (@$seq_ids) {
    my $SeqObj = &seq_id_to_SeqObj($dbh, $seqid);
    &add_features_to_SeqObj($dbh, $seqid, $SeqObj);
    $OUT->write_seq($SeqObj);
}
