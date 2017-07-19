#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

my $opt = {};
&getopts('D:i:n:o:s:S', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};;
my $seqid = $opt->{s};
my $split = $opt->{S};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'access');

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $OUT;
if ($outfile && !$split) {
    my $o = $outfile;
    $OUT = Bio::SeqIO->new(-file => ">$o",
			   -format => 'genbank');
}
if (! $outfile) {
    $OUT = Bio::SeqIO->new(-fh => \*STDOUT,
			   -format => 'genbank');
}



my $seq_ids = $seqid ? [$seqid] : &set_id_to_seq_ids($dbh, $setid);
foreach my $seqid (@$seq_ids) {
    if ($split) {
	my $o = $outfile . ".$seqid";;
	$OUT = Bio::SeqIO->new(-file => ">$o",
			       -format => 'genbank');
    }
    my $SeqObj = &seq_id_to_SeqObj($dbh, $seqid);
    &add_features_to_SeqObj($dbh, $seqid, $SeqObj);
    $OUT->write_seq($SeqObj);
}
