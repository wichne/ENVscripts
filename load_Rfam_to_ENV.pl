#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use Bio::SeqFeature::Generic;
use DBI;

my $args = {};
&getopts('i:d:u:p:', $args);

my $input_file=$args->{i};
my $db = $args->{d};
my $pwd = $args->{p};
my $user = $args->{u} ? $args->{u} : $ENV{USER};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $user, $pwd);

open (my $IN, $input_file) or die "Can't open $input_file: $!\n";

my %RFAM;
my %SEQ;
while (my $line = <$IN>) {
    next if ($line =~ /^#/);
    chomp $line;
    # Okay, somewhere the rfam_scan output changed from a gff to a columnar output, so new parsing 
#    my ($seq_acc,
#	$source,
#	$method,
#	$start,
#	$end,
#	$score,
#	$dir,
#	$frame,
#	$attributes) = split /\s+/, $line;
    my ($seq_acc, $blank_accession,
	$rfam_name, $rfam_acc, $model_type, $mdl_from, $mdl_to,
	$seq_from, $seq_to, $dir,
	$truncated, $pass, $gc, $bias,
	$bit_score, $evalue, $inc, $description) = split /\s+/, $line;
    my ($min, $max) = sort {$a<=>$b} ($seq_from, $seq_to);
    my $strand = $dir eq "-" ? -1 : 1;
    my ($type, $name, $gene, $so_term) = &annotation_from_Rfam($dbh, $rfam_acc);

    my @accs = split(/\|/, $seq_acc);
    my ($seq_id, $seqobj, $sacc);
    foreach my $ax(@accs) {
	if (!$SEQ{$ax}) {
#	    $seq_id = get_seq_id_by_seq_accession($dbh, $ax);
	    $seqobj = get_sequence_by_accession($dbh, $ax);
	    $sacc = $ax;
	    if ($seqobj) { 
		$SEQ{$ax} = $seqobj;
		$seq_id = $seqobj->display_id;
		last;
	    }
	} else {
	    $seqobj = $SEQ{$ax};
	    $seq_id = $seqobj->display_id;
	    last;
	}
    }
    if (!$seq_id) { die " Couldn't get a sequence_id from the accession '$seq_acc'.\n" }

    # does this feature already exist?
    my $coords_q = "SELECT sfm.feature_id FROM seq_feat_mappings sfm, sequence_features f"
	. " WHERE sfm.seq_id = $seq_id"
	. " AND sfm.strand = \"$strand\""
	. " AND (($min >= sfm.feat_min AND $min <= sfm.feat_max)"
	. " OR ($max >= sfm.feat_min AND $max <= sfm.feat_max)"
	. " OR ($min <= sfm.feat_min AND $max >= sfm.feat_max))"
	. " AND f.feature_id=sfm.feature_id and f.feat_type = \"$type\"";
    my $feat_r = $dbh->selectcol_arrayref($coords_q);
    if (@$feat_r) {
	foreach my $fid (@$feat_r) {
	    my $featr = get_features_by_feature_id($dbh, $fid);
	    if ($featr->{$fid}->{'location'}->{$seq_id}->{'feat_min'} != $min ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'feat_max'} != $max ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'strand'} != $strand #||
#		$featr->{$fid}->{'location'}->{$seq_id}->{'min_partial'} != $min_partial ||
#		$featr->{$fid}->{'location'}->{$seq_id}->{'max_partial'} != $max_partial
		) {
		update_feature_mapping($dbh, $seq_id, $fid, {'feat_min'    => $min,
							     'feat_max'    => $max,
							     'strand'      => $strand,
#							     'min_partial' => $min_partial,
#							     'max_partial' => $max_partial
				       });
	    }
	}
#	delete_feature_annotations($dbh, $fid);
#       update the annotation...
    } else {
	my $feature_acc = $sacc . "-" . $rfam_name . "-" . ++$RFAM{$rfam_name};
	
	my $featobj = new Bio::SeqFeature::Generic ( -start => $min,
						     -end => $max,
						     -strand => $strand,
						     -primary => $type,
						     -source_tag   => 'Rfam',
						     -display_name => $feature_acc,
						     -score  => $bit_score );
    
	my $seqobj;
	my $fid = load_SeqFeature($dbh, $seq_id, $featobj, $seqobj, $so_term);
	
	my $annObj = Bio::Annotation::Collection->new();
	$annObj->add_Annotation("product",
				Bio::Annotation::SimpleValue->new(-value=>$name));
	$annObj->add_Annotation("gene",
				Bio::Annotation::SimpleValue->new(-value=>$gene)) if ($gene);
	
	load_feature_annotations($dbh, $fid, $annObj, "Rfam", 3);
    }
}

sub annotation_from_Rfam {
    my $dbh = shift;
    my $model = shift;
    my $query = "SELECT feature_type, product, gene, SO_term from rfam.rfam where accession=\"$model\"";
    my $row = $dbh->selectrow_arrayref($query);
    if ($row) {
	return (@$row);
    } else {
	warn "No entry for $model in rfam db.\n";
    }
}
