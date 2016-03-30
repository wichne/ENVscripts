#!/usr/bin/perl
#Edited by JMo 11/07/2015. Delete old KO features before update
#Script takes 2-column output from Koala ko assignment of proteins from genomes 
#and inserts into the evidence_table.


use Getopt::Std;
use DBI;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

my %arg;
&getopts('D:u:p:i:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
if (! $infile) { die "No KOALA file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}

open my $in, $infile or die "Can't open $infile: $!\n";

while (my $line = <$in>) {
    chomp $line;
    my ($acc, $ko) = split/\s+/, $line;
    next if (! $ko);
    my @acc = $acc;
    my $fid_ref = &get_feature_id_by_accession($dbh, $acc);
    my $fid = $fid_ref->{$acc}->{'feature_id'};

    my $fref = &get_protein_by_feature_id($dbh, $fid);
    my $ev_d = "DELETE FROM feature_evidence WHERE feature_id=$fid AND ev_type='KO'";
    my $dsth= $dbh->do($ev_d);
    my $ev_i = "INSERT INTO feature_evidence (feature_id, ev_accession, feat_min, feat_max, ev_type, program)"
	. " VALUES ($fid, \"$ko\", 1, " . length($fref->{$fid}->{'product'}) . ", \"KO\", \"KOALA\")\n";

    my $sth = $dbh->do($ev_i);
}

##subs##
sub get_feature_id_by_accession {
    my $dbh = shift;
    my @accession = @_;

    my $f_id_q = "SELECT accession, feature_id FROM feature_accessions"
        . " WHERE accession in (\"" . join("\",\"", @accession) . "\")";

    my $r = $dbh->selectall_hashref($f_id_q, 'accession'); 
    {
        return $r;
    }
}

sub get_protein_by_feature_id {
    my $dbh = shift;
    my @feature_ids = @_;

    # grab all the sequences based on feature_id
    my $prot_q = "SELECT sf.feature_id, sf.product FROM sequence_features sf JOIN seq_feat_mappings sm using (feature_id)"
        . " WHERE feature_id in (" . join(",",@feature_ids) . ")";

    my $sequences = $dbh->selectall_hashref($prot_q, "feature_id");

    return $sequences;
}
