#!/usr/bin/perl

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
if (! $infile) { die "No input file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}

open my $in, $infile or die "Can't open $infile: $!\n";

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " SELECT feature_id , ?, ?, 'Interpro', ?, ?, ?, ?, ?, ?"
    . " FROM feature_accessions WHERE accession=?";
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    chomp $line;
    my @f = split/\t/, $line;
    if (@f < 11) { die "Input file is not in Interpro tsv format\n"; }
    if ($f[9] ne "T") { next }
    my $acc;
    if ($f[0] =~ /\|/) {
	(undef, $acc, undef) = split/\|/,$f[0],3;
    } else {
	$acc = $f[0];
    }
    
    $sth->execute($f[6], $f[7], $f[3], $f[4], undef, undef, undef, $f[8], $acc);

    if ($f[11]) {
	$sth->execute($f[6], $f[7], 'iprscan', $f[11], undef, undef, undef, $f[8], $acc);
    }
}
