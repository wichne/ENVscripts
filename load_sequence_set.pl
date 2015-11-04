#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;

&getopts('D:u:p:i:s:');
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", $opt_u, $opt_p);

my $setname = $opt_s;

my $set_id = &load_sequence_sets($dbh, $setname, undef, undef);

open my $in, $opt_i;

my @seq_ids;
while (my $l = <$in>) {
    chomp $l;
    push @seq_ids, &get_seq_id_by_accession($dbh, $l);
}
&link_seq_to_set($dbh, $set_id, @seq_ids);
