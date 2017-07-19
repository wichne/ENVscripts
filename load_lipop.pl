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

my $ev_d = "DELETE FROM feature_evidence"
    . " WHERE feature_id = ?"
    . " AND ev_type='SP'"
    . " AND program='LipoP-1.0'";
my $dsth = $dbh->prepare($ev_d);

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, score)"
    . " VALUES(?, 1, ?, 'LipoP-1.0', 'SP', 'SpII', ?)";
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    next if ($line =~ /^\#/);
    chomp $line;
    $line =~ s/^\s*//;
    my ($acc,
	$ev_type,
	$score,
	$margin,
	$cleavage,
	$plus2) = split/\s+/, $line;
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

     if ($ev_type eq "SpII"){
	my $str = "$score;$margin;$cleavage;$plus2;";
	$cleavage =~ /(\d+)/;
	my $cpos = $1;
	my @acc= split(/\|/, $acc);
	$sth->execute($fid, $cpos, $str);
    }
}
