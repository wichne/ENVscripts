#!/usr/bin/perl
#edited by JMo July 16 2014
use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use DBI;

my $host = $ENV{'DBSERVER'} ? $ENV{'DBSERVER'} : 'localhost';

my %arg;
&getopts('D:P:U:i:h',\%arg);
my $USER = $arg{'U'} ? $arg{'U'} : $ENV{'USER'}; 
my $PSWD = $arg{'P'} or die "Need to provide database password with -P\n";
my $DB = $arg{'D'} or die "Need to provide database with -D\n";

if ($arg{'h'}) {
    print "USAGE load_dbCAN_to_ENV.pl -D db -P dbpassword -i inputfile [-U user ]\n";
}

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$DB", $USER, $PSWD);

my $infile = $arg{'i'} or die "Need to provide input file with -i\n";

open (my $IN, $infile);

my $delete_q = "DELETE e FROM feature_evidence e, feature_accessions a"
    . " WHERE a.accession = ?"
    . " AND e.feature_id=a.feature_id"
    . " AND e.accession = ?";
my $insert_q = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " SELECT a.feature_id, ?, ?, \"dbCAN\", \"CAZy\","
    . " ?, ?, ?, ?, ?"
    . " FROM feature_accessions a"
    . " WHERE a.accession = ?";
my $sth = $dbh->prepare($insert_q);
my $linecount = 0;
while (my $line = <$IN>) {
    if ($line =~ /^#/) { next }
    $linecount++;
    chomp $line;
    my @f = split/\s+/, $line;
#    my @acc = split/\|/, $f[3]; #hmmscan of pfam in domtblout with multiple acc, be sure to change subarray in execute if needed
#    my $acc = $f[4] eq "-" ? $f[3] : $f[4]; #use this line if hmmscan and have single acc in field 3 or 4
    my @acc = split/\|/, $f[0]; #use this line if hmmsearch and acc in field 1
    my $tacc;
    foreach my $x (@acc) {
	if ($x ne "gnl" && $x ne "fig") {
	    $tacc = $x;
	    last;
	}
    }
    if (! $tacc) { warn "Can't parse protein accession from '$line'. Not loading this one.\n"; next } 

    my $hacc = $f[4];
    $hacc =~ s/\.\d+$//;
    my $sth = $dbh->prepare($insert_q);
    my $sxs = $sth->execute($f[17], $f[18], $hacc, $f[15], $f[16], $f[5],
			    $f[13], $tacc);
    if (! $sxs) { print STDERR $dbh->errstr() }
}
