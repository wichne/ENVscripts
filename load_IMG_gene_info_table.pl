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
    print STDERR "load_IMG_gene_info_table.pl -i [inputfile] -D [ENVdatabase] -p [dbpswd]

This program loads particular information from the tab-delimited table you get when you export
gene features from your Gene Cart into the feature_annotations and feature_evidence tables in 
an ENV-schema database.

Datatypes loaded: Product_name, gene_sym, EC

The program does NOT delete or update any existing data in the database. It only inserts new rows.
It does, however, avoid inserting duplicate KO accessions.

Accessions that are missing or are duplicate in the feature_accessions table are reported in
separate logfiles (load_IMG_gene_info_table.[pid].missing and load_IMG_gene_info_table.[pid].duplicate

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
	$gene_sym,
	$product,
	$seq_len,
	$prot_len,
	$organism,
	$ec, # prefixed with EC:, 
	$COGs, 
	$COG_cats, # letters
	$Pfams, # space separated
	$TIGRfams,
	$sigp, # Yes or No?
	$tmh, # number or No
	$genome_id) = split/\t/, $line;

    if ($oid eq "gene_oid") { next } # this is the header line

    $DATA{$locus}->{'product'} = $product;
    while ($ec =~ /^EC\:(\S+)/g) {
	push @{$DATA{$locus}->{'ec'}}, $1;
    }
    $DATA{$locus}->{'gene_sym'} = $gene_sym if ($gene_sym);
#    while ($COGs =~ /^(COG\d{4})/g) {
#	push @{$DATA{$locus}->{'COG'}}, $1;	
#    }
    $DATA{$locus}->{'protlen'} = $prot_len;
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
    open my $out, ">load_IMG_gene_info_table.$$.missing";
    print $out $missing;
    close $out;
    warn "$mc accessions were missing from the db - see load_IMG_gene_info_table.$$.missing for details.\n";
}
if ($duplicate) {
    open my $out, ">load_IMG_gene_info_table.$$.duplicate";
    print $out $duplicate;
    close $out;
    warn "$dc accessions are assigned to more than one feature_id in the db - see load_IMG_info.$$.duplicate for details.\n";
}

# Now load the info
my $anni = "INSERT INTO feature_annotations"
    . " (feature_id, data_type_id, value, rank, source, edit_by)"
    . " VALUES (?, ?, ?, '$rank', '$source', USER())";
my $annsth = $dbh->prepare($anni);

#my $evi = "INSERT INTO feature_evidence"
#    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
#    . " VALUES (?, ?, ?, 'IMG_pipeline', ?, ?, ?, ?, ?, ?)";
#my $evsth = $dbh->prepare($evi);

foreach my $acc (keys %DATA) {
    my $fid = $DATA{$acc}->{'feature_id'};
    if (!$fid) { next }
    $annsth->execute($fid, 66, $DATA{$acc}->{'product'});
    $annsth->execute($fid, 35, $DATA{$acc}->{'gene_sym'}) if (defined $DATA{$acc}->{gene_sym});
    $annsth->execute($fid, 1, join(" ", @{$DATA{$acc}->{'ec'}})) if (defined $DATA{$acc}->{'ec'});
    
#    if (defined $DATA{$acc}->{'COG'}) {
#    foreach my $row (@{$DATA{$acc}->{'COG'}}) {
#	my ($evacc, $evalue) = @$row;
#	$evsth->execute($fid, 1, $DATA{$acc}->{'protlen'}, 'COG', $evacc, 0, 0, 0, $evalue);
#    }
}

#    if (defined $DATA{$acc}->{'KO'}) {
#    foreach my $row (@{$DATA{$acc}->{'KO'}}) {
#	my ($evacc, $evalue) = @$row;
#	if (&check_ko($dbh, $fid, $evacc)) {
#	    $evsth->execute($fid, 1, $DATA{$acc}->{'protlen'}, 'KO', $evacc, 0, 0, 0, $evalue);
#	}
#    }}
#}

sub check_ko {
    my $dbh = shift;
    my $fid = shift;
    my $koacc = shift;

    my $q = "SELECT count(*) from feature_evidence where feature_id=$fid AND ev_accession=\"$koacc\"";
    my @r = $dbh->selectrow_array($q);
    if ($r[0] > 0) { return 0 } else { return 1 }
}
