#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;
use DBI;

my %arg;
&getopts('i:D:p:h', \%arg);
my $source = "IMG";
my $rank = 9;

if ($arg{'h'}) {
    print STDERR "load_IMG_info.pl -i [inputfile] -D [ENVdatabase] -p [dbpswd]

This program loads particular information from a tab-delimited version of the IMG .info file
into the feature_annotations and feature_evidence tables in an ENV-schema database.

Datatypes loaded: Product_name, EC, COG, KO

The program does NOT delete or update any existing data in the database. It only inserts new rows.
It does, however, avoid inserting duplicate KO accessions.

Accessions that are missing or are duplicate in the feature_accessions table are reported in
separate logfiles (load_IMG_info.[pid].missing and load_IMG_info.[pid].duplicate

";
    exit;
}

my $file = $arg{'i'};
my $db = $arg{'D'};
my $pswd = $arg{'p'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $ENV{USER}, $pswd);

open(my $in, $file) or die "Can't open $file for read: $!\n";

# Start by reading data into a structure
my %DATA;
while (my $line = <$in>) {
    chomp $line;
    my ($oid,
	$locus,
	$source,
	$cluster_info,
	$gene_info,
	$evalue) = split/\t/, $line;
    if ($oid eq "gene_oid") { next }
    if ($source eq "Product_name") {
	$DATA{$locus}->{'product'} = $gene_info;
    } elsif ($source =~ /^EC\:(\S+)/) {
	push @{$DATA{$locus}->{'ec'}}, $1;
    } elsif ($source =~ /^(COG\d{4})/) {
	push @{$DATA{$locus}->{'COG'}}, [$1, $evalue];	
    } elsif ($source =~ /^KO\:(K\d{5})/) {
	push @{$DATA{$locus}->{'KO'}}, [$1, $evalue];
    } elsif ($source eq "Protein_length") {
	$gene_info =~ s/aa//;
	$DATA{$locus}->{'protlen'} = $gene_info;
    }
}

# Test for missing data
my $fidq = "SELECT feature_id FROM feature_accessions WHERE accession=?";
my $fsth = $dbh->prepare($fidq);
my $missing; my $mc;
my $duplicate; my $dc;
foreach my $acc (keys %DATA) {
    $fsth->execute($acc);
    my $rows = $fsth->fetchall_arrayref;
    if (scalar @$rows == 0) { 
	$missing .= "$acc\n";
	$mc++;
    } elsif (scalar @$rows > 1) {
	$duplicate .= "$acc";
	foreach my $idr (@$rows) {
	    $duplicate .= "\t$idr->[0]";
	}
	$duplicate .= "\n";
	$dc++;
    } else {
	$DATA{$acc}->{'feature_id'} = $rows->[0]->[0];
    }
}
if ($missing) {
    open my $out, ">load_IMG_info.$$.missing";
    print $out $missing;
    close $out;
    warn "$mc accessions were missing from the db - see load_IMG_info.$$.missing for details.\n";
}
if ($duplicate) {
    open my $out, ">load_IMG_info.$$.duplicate";
    print $out $duplicate;
    close $out;
    warn "$dc accessions are assigned to more than one feature_id in the db - see load_IMG_info.$$.duplicate for details.\n";
}

# Now load the info
my $anni = "INSERT INTO feature_annotations"
    . " (feature_id, data_type_id, value, rank, source, edit_by)"
    . " VALUES (?, ?, ?, '$rank', '$source', USER())";
my $annsth = $dbh->prepare($anni);

my $evi = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " VALUES (?, ?, ?, 'IMG_pipeline', ?, ?, ?, ?, ?, ?)";
my $evsth = $dbh->prepare($evi);

foreach my $acc (keys %DATA) {
    my $fid = $DATA{$acc}->{'feature_id'};
    if (!$fid) { next }
    $annsth->execute($fid, 66, $DATA{$acc}->{'product'});
    $annsth->execute($fid, 1, join(" ", @{$DATA{$acc}->{'ec'}})) if (defined $DATA{$acc}->{'ec'});
    
    if (defined $DATA{$acc}->{'COG'}) {
    foreach my $row (@{$DATA{$acc}->{'COG'}}) {
	my ($evacc, $evalue) = @$row;
	$evsth->execute($fid, 1, $DATA{$acc}->{'protlen'}, 'COG', $evacc, 0, 0, 0, $evalue);
    }}

    if (defined $DATA{$acc}->{'KO'}) {
    foreach my $row (@{$DATA{$acc}->{'KO'}}) {
	my ($evacc, $evalue) = @$row;
	if (&check_ko($dbh, $fid, $evacc)) {
	    $evsth->execute($fid, 1, $DATA{$acc}->{'protlen'}, 'KO', $evacc, 0, 0, 0, $evalue);
	}
    }}
}

sub check_ko {
    my $dbh = shift;
    my $fid = shift;
    my $koacc = shift;

    my $q = "SELECT count(*) from feature_evidence where feature_id=$fid AND ev_accession=\"$koacc\"";
    my @r = $dbh->selectrow_array($q);
    if ($r[0] > 0) { return 0 } else { return 1 }
}
