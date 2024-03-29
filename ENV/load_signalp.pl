#!/usr/bin/perl

use Getopt::Std;
use DBI;
use strict;
use lib $ENV{ENVSCRIPTS};
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

my $ev_d = "DELETE FROM feature_evidence"
    . " WHERE feature_id = ?"
    . " AND ev_type='SP'"
    . " AND program='signalp-4.1'";
my $dsth = $dbh->prepare($ev_d);

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, score)"
    . " VALUES (?, 1, ?, 'signalp-4.1', 'SP', 'SpI', ?)"
    ;
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    if ($line =~ /^\#/) { next } # skip header lines
    chomp $line;
    my @f = split/\s+/, $line;
    if (@f != 12) { die "Input file $infile is not in signalp-4.1 summary format: ". @f . "\n"; }

    my ($acc,
	$cmax,
	$cpos,
	$ymax,
	$ypos,
	$smax,
	$spos,
	$smean,
	$d,
	$sp,
	$dmaxcut,
	$networks)= @f;

    my $fid;
    if ($acc =~ /^\d+$/) {
	$fid = $acc;
    } else {
	$fid = get_feature_id_by_accession($dbh,$acc);
	if (!$fid) {
	    print "WARNING: Couldn't get a feature_id from accession '$acc'. Skipping.\n";
	    next;
	}
    }

    $dsth->execute($fid);
    if ($sp eq "Y") {
	my $score = "Cmax=$cmax;Ymax=$ymax;Smean=$smean;D=$d;Dmaxcut=$dmaxcut;";
	my @acc = split(/\|/, $acc);
	$sth->execute($fid, $ypos, $score);
    }
}
