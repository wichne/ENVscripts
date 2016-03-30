#!/usr/bin/perl
use DBI;
use lib "/share/scripts";
use ENV;
use Getopt::Std;
use strict;

my %arg;
&getopts('D:u:p:',\%arg);
my $db = $arg{D};
my $dbp= $arg{p};
my $dbu = $arg{u} ? $arg{u} : $ENV{USER};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $dbu, $dbp);

my %DATA;
my %COF;
my $cofactor_q = "SELECT s.name, cofactor, count(distinct e.feature_id)"
#    . " FROM feature_annotations a, seq_feat_mappings m, sequence_sets s, seq_set_link l, egad.ec_cofactor_link e"
    . " FROM feature_evidence e, seq_feat_mappings m, sequence_sets s, seq_set_link l, egad.hmm_cofactor_link h"
    . " WHERE s.set_id=l.set_id and m.seq_id=l.seq_id"
    . " AND e.feature_id=m.feature_id"
    . " AND (e.ev_accession=h.hmm_acc OR LEFT(e.ev_accession,7)=h.hmm_acc)"
    . " GROUP by s.name, cofactor";
my $sth = $dbh->prepare($cofactor_q);
$sth->execute();
while (my @row = $sth->fetchrow_array) {
    $DATA{$row[0]}->{$row[1]} = $row[2];
    $COF{$row[1]} = 1;
}

print "Genome\t";
print join("\t", sort keys %COF),"\n";

foreach my $genome (sort keys %DATA) {
    print "$genome";
    foreach my $cof (sort keys %COF) {
	print "\t$DATA{$genome}->{$cof}";
    }
    print "\n";
}
exit;
