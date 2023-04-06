#!/usr/bin/perl
# edited by JMo July 16 2014
# edited by WCN Feb 10, 2015
use lib $ENV{ENVSCRIPTS};
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
    print "USAGE load_hmms_to_ENV.pl -D db -P dbpassword [-U user ]\n";
}

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$DB", $USER, $PSWD);

my $infile = $arg{'i'} or die "Need to provide domtblout file with -i\n";

open (my $IN, $infile);

my $delete_q = "DELETE e FROM feature_evidence e"
    . " WHERE feature_id = ?"
    . " AND ev_type = 'HMM'"
    . " AND e.ev_accession = ?";
my $insert_q = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " VALUES(?, ?, ?, \"hmmerv3\", \"HMM\","
    . " ?, ?, ?, ?, ?)";
my $sth = $dbh->prepare($insert_q);
my $linecount = 0;
my $type;
while (my $line = <$IN>) {
    if ($line =~ /^#/) {
	# determine whether we're looking at hmmscan or hmmsearch output by examining header
	if ($line =~ /^# target/) {
	    $type = "scan";
	} elsif ($line =~ /^# query/) {
	    $type = "search";
	}
	next;
    }
    $linecount++;
    chomp $line;
    my @f = split/\s+/, $line;

    my (@acc,$tacc,$hacc);
    my ($acc_field, $hmm_field);
    if ($f[0] =~ /(TIGR|PF)\d{5}/) {
	$hacc = $f[0];
	$acc_field = $f[3];
    } elsif ($f[1] =~ /(TIGR|PF)\d{5}/) {
	$hacc = $f[1];
	$acc_field = $f[3];
    } elsif ($f[4] =~ /(TIGR|PF)\d{5}/) {
	$hacc = $f[4];
	$acc_field = $f[0];
    } else {
	die "Can't find HMM acc in $line\n";
    }


    @acc = split/\|/,$acc_field;

    my $tacc;
    foreach my $x (@acc) {
	if ($x eq "gnl" || $x eq "fig" || $x eq "gp" || $x eq $DB) { next }
	else { $tacc = $x; last; }
    }

    my $fid;
    if ($tacc =~ /^\d+$/) {
	$fid = $tacc;
    } else {
	$fid = get_feature_id_by_accession($dbh, $tacc);
	if (!$fid) {
	    print "WARNING: Couldn't get a feature_id from accession '$tacc'. Skipping.\n";
	    next;
	}
    }

    $hacc =~ s/\.\d+$//;

    $dbh->do($delete_q, {}, ($fid, $hacc));

    my $sth = $dbh->prepare($insert_q);
    my $sxs = $sth->execute($fid, $f[17], $f[18], $hacc, $f[15], $f[16], $f[5],
			    $f[13]);
    if (! $sxs) { print STDERR $dbh->errstr() }
}
