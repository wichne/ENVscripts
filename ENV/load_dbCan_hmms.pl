#!/usr/bin/perl

#09/16/2014, Jennifer Mobberley (mobberley.jennifer@pnnl.gov)
#This script will load the best-hit dbCAN hmm generated from local annotation scripts (see README file at csbl.bmb.uga.edu/dbCAN/download). An additional step to append dbCAN hmm model length was carried out using the /home/mobb021/hmm/dbCAN/dbCAN_hmm_prep.sh

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
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length,score)"
    . " SELECT feature_id , ?, ?, 'dbCAN', 'CAZy', ?, ?, ?, ?, ?"
    . " FROM feature_accessions WHERE accession=?";
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    chomp $line;
    my @f = split/\s+/, $line;
    if (@f != 9) { die "Input file $infile is not in dbcan summary format: ". @f . "\n"; }

    my ($acc,
	$ev_accession,
	$evalue,
	$ev_min,
	$ev_max,
	$feat_min,
	$feat_max,
	$cover,
	$ev_length)= @f;
    my $score = "evalue=$evalue;coverage=$cover;";
	$sth->execute($feat_min, $feat_max, $ev_accession, $ev_min, $ev_max, $ev_length, $score, $acc);
    
}
