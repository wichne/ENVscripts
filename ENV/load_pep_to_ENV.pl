#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use strict;
use lib $ENV{SCRIPTS};
use Getopt::Std;
use ENV;

my %arg;
&getopts('D:u:p:f:s:',\%arg);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg{D}", $arg{u}, $arg{p});

my $fo = Bio::SeqIO->new(-file => $arg{f},
			 -format=> 'fasta') || die "Bad fasta file: $arg{f}\n";
my $source = $arg{s} ? $arg{s} : "unspecified";

while (my $reco = $fo->next_seq) {
    my $id = $reco->display_id;
    my $pep = $reco->seq;
    $pep =~ s/\*$//;

    my @acc = split /\|/,$id;
    my $accs = "'" . join("','", @acc) . "'";
    my $acc_q = "SELECT distinct a.feature_id, length(product) FROM feature_accessions a, sequence_features f"
	. " WHERE accession in ($accs)"
	. " AND f.feature_id=a.feature_id";
    my $r = $dbh->selectall_arrayref($acc_q);
    if (@$r == 0) {
	# First, see if underlying contig is loaded
	# If it's not, throw a warning and move on.
	$id =~ /(.*)\_\d+$/;
	my $seqacc = $1;
	my $seq_id = ENV::get_seq_id_by_accession($dbh,$seqacc);
	if (! $seq_id) { 
	    warn "No seq_id for $seqacc (from $id)";
	    next;
	} else {
	    warn "Loading $id onto $seqacc ($seq_id)";
	}

	# parse out needed information from description line:
	my (undef, $min, $max, $strand, $desc) = split / *# +/, $reco->description;
	my @ann = split/;/,$desc;

	# load the sequence_feature record
	my $sfi = "INSERT INTO sequence_features (feat_type, product, inserted_by, date_inserted) VALUES ('CDS', \"$pep\", USER(), CURDATE())";
	$dbh->do($sfi);
	my $fid = $dbh->last_insert_id("%", "%", "", "");

	# Now load the seq_feat_mappings record
	my $sfmi = "INSERT INTO seq_feat_mappings"
	    . " (seq_id, feature_id, feat_min, feat_max, strand, phase)"
	    . " VALUES ($seq_id, $fid, $min, $max, \"$strand\", 0)";
	$dbh->do($sfmi);
	
	# load the accession
	my $fai = "INSERT INTO feature_accessions"
	    . " (feature_id, source, accession)"
	    . " VALUES ($fid, '$source', '$id')";
	$dbh->do($fai);

	# and here we would load annotation

    } elsif (@$r == 1) {
	if ($r->[0]->[1] == 0) {
	    my $q = "UPDATE sequence_features"
		. " SET product=\"$pep\""
		. " WHERE feature_id=$r->[0]->[0]";
	    $dbh->do($q);
	}
    } else {
	foreach my $s(@$r) {
	    print "!!$s->[0]\t$s->[1]\n";
	}
    }
}
