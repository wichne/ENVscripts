#!/usr/bin/perl
#Updated JMo 06/13/2015, modified acc to accomodate new peptide files

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
    . " AND ev_type='TMH'"
    . " AND program='TMHMM2.0'";
my $dsth = $dbh->prepare($ev_d);

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession)"
    . " VALUES(?, ?, ?, 'TMHMM2.0', 'TMH', 'TMhelix')";
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    next if ($line =~ /^\#/);
    chomp $line;
    my $short = 0;
    if ($line =~ /ExpAA/) {
	$short = 1;
    }
    my @f = split/\t/, $line;
#    if (@f != 6) { die "Input file is not in TMHMM2.0 summary format\n"; }
    
    my $acc = $f[0];
    my $fid;
    if ($acc =~ /^\d+$/) {
	$fid = $acc;
    } else {
	$fid = get_feature_id_by_accession($dbh, $acc);
	if (!$fid) {
	    print "WARNING: Couldn't get a feature_id from accession '$acc'. Skipping.\n";
	    next;
	}
    }

    $dsth->execute($fid);

    if ($short) {
	my (undef, $len) = split/\=/, $f[1];
	my (undef, $string) = split/\=/, $f[5];
	
	while ($string =~ /(\d+)\-(\d+)/g) {
	    $sth->execute($fid, $1, $2);
	}
    } else {
	# this is for long form output.
	if ($f[2] eq "TMhelix") {
	    #my @acc = split(/\|/, $acc);
	    my $ev_acc = $f[0];
	    my $min = $f[3];
	    my $max = $f[4];
	    $sth->execute($fid, $min, $min, $ev_acc);
	}
    }
}
