#!/usr/bin/perl -w

use DBI;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;

# input args
my $args = {};
&getopts('f:F:D:u:p:h', $args);
my $filename = $args->{'f'} or die "Need to provide filename with -f\n";
my $format = $args->{'F'} or die "Provide file format (BBH, OrthoFinder) with -F\n";

# connect to db
my $dbh;
my ($user, $db, $password);
if (defined $args->{'u'} && $args->{'u'}) { $user = $args->{'u'} }
else { $user = $ENV{USER}; }
if (defined $args->{'p'}) { $password = $args->{'p'} }
if (defined $args->{'D'}) { $db = $args->{'D'} }
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
$dbh = connect($args);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}
#print STDERR "Connected to $host:$db as $user\n";

# open input file (build_ortholog_table output)
# a tab-delimited file
open (my $inh, $filename);

my $orthfam_i = "INSERT INTO ortho_group (feature_id, ortho_group, method)"
    . " SELECT distinct feature_id, ?, ?"
    . " FROM feature_accessions a"
    . " WHERE (accession = ? OR feature_id = ?)";
my $orthfam_h = $dbh->prepare($orthfam_i);

if ($format eq 'BBH') {
    # first line of file should be a header, and the genome fields should match the
    # common_name field in the 'organisms' table in the db.
    my $line = <$inh>;
    chomp $line;
    # do a quick data check to see if this file is in the right format
    if ($line !~ /^MCLFamID\b/) { die "First line indicates this file is not in the correct format for this program (or header line is missing).\n"; }
    my @fields = split(/\t/, $line);
    my @genomes = @fields[11..$#fields];

    # now read in each line and process it
    while ($line = <$inh>) {
	chomp $line;
	my ($mclid,
	    $bbhid,
	    $synfamid,
	    $tigrfamid,
	    $genomes,
	    $genes,
	    $product,
	    $gene,
	    $roleid,
	    $ec,
	    $tc,
	    @dbs) = split(/\t/, $line);

	for (my $dbidx=0; $dbidx<@genomes;$dbidx++) {
	    if (! $dbs[$dbidx]) { next }
	    my @genes = split(/\s+/,$dbs[$dbidx]);
	    foreach my $g (@genes) {
		if ($mclid) {
		    my $mclok = $orthfam_h->execute($mclid, "MCL", $g, $g);
		    if (! defined $mclok) { die "Darn. " . $orthfam_h->{Statement}, "\n"; }
		}
		if ($bbhid) {
		    my $bbhok = $orthfam_h->execute($bbhid, "BBH", $g, $g);
		    #if (! defined $bbhok) { warn "Darn. " . $orthfam_h->{Statement}, "\n"; }
		}
		if ($synfamid) {
		    my $synok = $orthfam_h->execute($synfamid, "synteny", $g, $g);
		    #if (! defined $synok) { die "Darn. " . $orthfam_h->{Statement}, "\n"; }
		}
		if ($tigrfamid) {
		    my $tigrok = $orthfam_h->execute($tigrfamid, "TIGRfam", $g, $g);
		    #if (! defined $tigrok) { die "Darn. " . $orthfam_h->{Statement}, "\n"; }
		}
	    }
	}
    }
} elsif ($format eq "OrthoFinder") {
    # first line of file should be a header
    my $line = <$inh>;
    chomp $line;
    # do a quick data check to see if this file is in the right format
    if ($line !~ /^Orthogroup\b/) { die "First line indicates this file is not in the $format format (or header line is missing).\n"; }
    my @fields = split(/\t/, $line);
    my @genomes = @fields[1..$#fields];

    while ($line = <$inh>) {
	chomp $line;
	my $orthid;
	my @orgProts;
	($orthid, @orgProts) = split(/\t/, $line);
	foreach my $orgList (@orgProts) {
	    my @prots = split(/, /, $orgList);
	    foreach my $id(@prots) {
		if ($id =~ /^((sp|tr)\|)?([^\|]+)\|*/) {
		    my $acc = $3;
		    #my $ofok = $orthfam_h->execute($orthid, "OrthoFinder", $acc);
		    #if (! defined $ofok) { warn "Darn $orthid\t$acc\n"; }
		    print "Inserting $orthid, $acc\n";
		} else {
		    die "$line\nI don't recognize this id '" . $id . "' as valid. Please update the program.\n";
		}
	    }
	}
    }
}
