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
    . " AND ev_type='PSORT'"
    . " AND program='PSORTb_v3.0'";
my $dsth = $dbh->prepare($ev_d);

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, score)"
    . " VALUES(?, 1, 1, 'PSORTb_v3.0', 'PSORT', ?, ?)";
my $sth = $dbh->prepare($ev_i);

while (my $line = <$in>) {
    next if ($line =~ /^\#/);
    chomp $line;
    my @f = split(/\t/, $line);
    if (@f == 3) { # terse output
	my ($acc, $loc, $score) = @f;
	$acc =~ s/\s+.*//; # field is fasta header, so there may be description present
	my $fid = get_feature_id_by_accession($dbh,$acc);
	if (!$fid) {
	    # see if this is a feature_id
	    if ($acc =~ /^\d+$/) {
		# put in a test here to see if it is a feature_id
		$fid = $acc;
	    }
	}
	if (!$fid) {
	    print "WARNING: Couldn't get a feature_id from accession '$acc'. Skipping.\n";
	    next;
	}
	
	# remove existing data
	$dsth->execute($fid);
	
	# insert result
	# print "For $acc ($fid) inserting '$loc', '$score'\n";
	$sth->execute($fid, $loc, $score);
    } elsif (@f == 35) { # long output
	my ($acc, $desc) = split(/\s/, $f[0], 2);
	my $loc = $f[30];
	$loc = $f[31] if ($loc eq "Unknown" and $f[31] =~ /\w+/);
	$loc .= "::" . $f[32] if ($f[32] =~ /\w+/);
	my $score = $f[32];
	if ($loc ne "Unknown") {
	    my $fid = get_feature_id_by_accession($dbh,$acc);
	    if (!$fid) {
		print "WARNING: Couldn't get a feature_id from accession '$acc'. Skipping.\n";
		next;
	    } else {
		$dsth->execute($fid);
		$sth->execute($fid, $loc, $score);
	    }
	}
    } else {
	warn "Wrong number of cols. Check file format?\n";
    }
}
    
