#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;
use DBI;

my %arg;
&getopts('i:D:p:h', \%arg);
my $source = "RAST";
my $rank = 9;

if ($arg{'h'}) {
    print STDERR "load_RAST_info.pl -i [inputfile] -D [ENVdatabase] -p [dbpswd]

This program loads particular information from the RAST 'Spreadsheet (tab-separated text format)'
into the feature_annotations and feature_evidence tables in an ENV-schema database.

Datatypes loaded: function (with EC parsed out), figfam

The program does NOT delete or update any existing data in the database. It only inserts new rows.

Accessions that are missing or are duplicate in the feature_accessions table are reported in
separate logfiles (load_RAST_info.[pid].missing and load_RAST_info.[pid].duplicate

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
    my ($ctg_id,
	$acc,
	$type,
	$location,
	$start,
	$stop,
	$strand,
	$function,
	$alias,
	$figfam,
	$evidence_codes,
	$nucseq,
	$protseq) = split/\t/, $line;
    if ($ctg_id eq "contig_id") { next }
    elsif (! &check_ctg_id($dbh, $ctg_id)) { die "contig $ctg_id is not in database: aborting.\n" }

    $type = "CDS" if ($type eq "peg");
    $acc =~ s/fig\|//;
    my @EC; my @TC;
    my $stripped_function = $function;
    while ($function =~ /(\s*\(EC (\d+(\.(\d+|-)){1,3})\))/g) {
	push @EC, $2;
	$stripped_function =~ s/\Q$1\E//;
    }
    while ($function =~ /(\s*\(TC (\d+\.[A-Z]\.\d+\.\d+\.\d+)\))/g) {
	push @TC, $2;
	$stripped_function =~ s/\Q$1\E//;
    }
    $DATA{$acc}->{'product'} = $stripped_function;
    $DATA{$acc}->{'ec'} = \@EC if (@EC);
    $DATA{$acc}->{'FIGFAM'} = $figfam if ($figfam);
    $DATA{$acc}->{'protlen'} = length($protseq);
    $DATA{$acc}->{'feature_id'} = &find_feature_id($dbh, $acc, $type, $ctg_id, $start, $stop, $strand, $protseq);
}

# Now load the info
my $anni = "INSERT INTO feature_annotations"
    . " (feature_id, data_type_id, value, rank, source, edit_by)"
    . " VALUES (?, ?, ?, '$rank', '$source', USER())";
my $annsth = $dbh->prepare($anni);

my $evi = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " VALUES (?, ?, ?, 'RAST_pipeline', ?, ?, ?, ?, ?, ?)";
my $evsth = $dbh->prepare($evi);

foreach my $acc (keys %DATA) {
    my $fid = $DATA{$acc}->{'feature_id'};
    if (!$fid) { next }
    $annsth->execute($fid, 66, $DATA{$acc}->{'product'});
    $annsth->execute($fid, 1, join(" ", @{$DATA{$acc}->{'ec'}})) if (defined $DATA{$acc}->{'ec'});
    
    if (defined $DATA{$acc}->{'FIGFAM'}) {
	$evsth->execute($fid, 1, $DATA{$acc}->{'protlen'}, 'FIGFAM',
			$DATA{$acc}->{'FIGFAM'}, 0, 0, 0, 1);
    }
}

sub find_feature_id {
    my $dbh = shift;
    my $acc = shift;
    my $type = shift;
    my $seq_acc = shift;
    my $end5 = shift;
    my $end3 = shift;
    my $strand = shift;
    my $protseq = shift;

    my ($low, $high) = $end5 < $end3 ? ($end5, $end3) : ($end3, $end5);

    # Test for missing data
    my $fid = &fid_from_acc($dbh, $acc);
    $fid = &fid_from_loc($dbh, $acc, $type, $seq_acc, $low, $high, $strand) if (! $fid);
    $fid = &insert_seqFeature($dbh, $acc, $seq_acc, $type, $low, $high, $strand, $protseq) if (! $fid);
    return $fid;
}

sub fid_from_acc {
    my $dbh = shift;
    my $acc = shift;
    my $fidq = "SELECT feature_id FROM feature_accessions WHERE source='RAST' AND accession=?";
    my $fsth = $dbh->prepare($fidq);
    $fsth->execute($acc);
    my $rows = $fsth->fetchall_arrayref;
    if (scalar @$rows > 1) {
	my $dups;
	foreach my $idr (@$rows) {
	    $dups .= "$idr->[0]::";
	}
	warn "$acc maps to multiple feature_ids: $dups\nUsing " . $rows->[0]->[0] . "\n";
    }
    print "$acc\t" . $rows->[0]->[0] . "\n" if (@$rows);
    return $rows->[0]->[0];
}

sub fid_from_loc {
    my $dbh = shift;
    my $acc = shift;
    my $type = shift;
    my $seq_acc = shift;
    my $low = shift;
    my $high = shift;
    my $dir = shift;
    my $strand = $dir eq "+" ? 1 : $dir eq "-" ? -1 : 0;
    my $fidq = "SELECT f.feature_id FROM seq_feat_mappings m, sequence_features f, sequence_accessions s"
	. " WHERE f.feat_type=\"$type\""
	. " AND s.seq_accession=?"
	. " AND m.seq_id=s.seq_id"
	. " AND m.feature_id=f.feature_id"
	. " AND m.strand = \"$strand\""
	. " AND (m.feat_min=$low OR m.feat_max=$high)";
    my $fsth = $dbh->prepare($fidq);
    $fsth->execute($seq_acc);
    my $rows = $fsth->fetchall_arrayref;
    if (scalar @$rows > 1) {
	my $dups;
	foreach my $idr (@$rows) {
	    $dups .= "$idr->[0]::";
	}
	warn "$acc maps to multiple feature_ids: $dups\nUsing " . $rows->[0]->[0] . "\n";
    }
    &insert_feature_accession($dbh, $acc, $rows->[0]->[0]) if (@$rows);
    print "$acc\t" . $rows->[0]->[0] . "\n" if (@$rows);
    return $rows->[0]->[0];
}

sub insert_seqFeature {
    my $dbh = shift;
    my ($acc, $seq_acc, $type, $low, $high, $strand, $protseq) = @_;
    my $dir = $strand eq "+" ? 1 : $strand eq "-" ? -1 : 0;
    my $phase = 0;
    my $min_partial = 0;
    my $max_partial = 0;
    my $sequence_features_i = "INSERT sequence_features"
	. " (feat_type_id, feat_type, product, inserted_by, date_inserted)"
	. " SELECT id, \"" . $type . "\", \"$protseq\", USER(), NOW()"
	. " FROM INSDC.feature_key"
	. " WHERE feature_key = \"" . $type . "\"";
#    print $sequence_features_i, "\n";
    my $row1 = $dbh->do($sequence_features_i);
    if (! defined $row1) {
	warn "\n!!!Couldn't execute '$sequence_features_i': $DBI::errstr\n";
	return;
    }
    
    my $feat_id = $dbh->last_insert_id("%", "%", "", "");
    my $seq_feat_mappings_i = sprintf "INSERT seq_feat_mappings"
	. " (seq_id, feature_id, feat_min, feat_max, strand,"
	. " phase, min_partial, max_partial)"
	. " SELECT seq_id, $feat_id, $low, $high, '$dir', '$phase', $min_partial, $max_partial"
	. " FROM sequence_accessions where seq_accession = \"$seq_acc\"";
#    print $seq_feat_mappings_i, "\n";
    my $row2 = $dbh->do($seq_feat_mappings_i);
    if (! defined $row2) {
	warn "\n!!!Couldn't execute '$seq_feat_mappings_i'\n";
	return;
    }
    &insert_feature_accession($dbh, $acc, $feat_id);
    print "$acc\t" . $feat_id . "\n";
    return $feat_id;
}

sub insert_feature_accession {
    my ($dbh, $acc, $feat_id) = @_;
    my $feat_acc_i = "INSERT feature_accessions"
	. " (feature_id, source, prefix, accession)"
	. " VALUES ($feat_id, \"$source\", \"fig\", \"$acc\")" ;
    my $row3 = $dbh->do($feat_acc_i);
 
}

sub check_ctg_id {
    my ($dbh, $ctg_id) = @_;
    my $q = "SELECT seq_id from sequence_accessions where seq_accession=\"$ctg_id\"";
    my $r = $dbh->selectall_arrayref($q);
    return scalar(@$r);
}
