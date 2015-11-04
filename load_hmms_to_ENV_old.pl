#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use DBI;

my $host = $ENV{'DBSERVER'} ? $ENV{'DBSERVER'} : 'localhost';

my %arg;
&getopts('D:P:U:i:',\%arg);
my $USER = $arg{'U'} ? $arg{'U'} : $ENV{'USER'}; 
my $PSWD = $arg{'P'} or die "Need to provide database password with -P\n";
my $DB = $arg{'D'} or die "Need to provide database with -D\n";

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$DB", $USER, $PSWD);

my $infile = $arg{'i'} or die "Need to provide input file with -i\n";

open (my $IN, $infile);

my $insert_q = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " SELECT a.feature_id, ?, ?, \"hmmerv3\", \"HMM\","
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
    my $acc = $f[4] eq "-" ? $f[3] : $f[4];
    my @tacc = split/\|/, $f[0];
    my $sth = $dbh->prepare($insert_q);
    my $sxs = $sth->execute($f[17], $f[18], $f[1], $f[15], $f[16], $f[5],
			    $f[13], $acc);
    if (! $sxs) { print STDERR $dbh->errstr() }
}
