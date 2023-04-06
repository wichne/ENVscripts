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

my %AA;
my %SEQ;
while (my $line = <$IN>) {
    chomp $line;
    my ($seq_acc,
	$index,
	$seq_from,
	$seq_to,
	$aa,
	$acodon,
	$intron_end5,
	$intron_end3,
	$score) = split /\s+/, $line;
    if ($index !~ /\d+/) { next }
    $AA{$aa}++;
    my ($min, $max, $strand) = $seq_from < $seq_to ? ($seq_from, $seq_to, 1) : ($seq_to, $seq_from, -1);
    my @accs = split/\|/, $seq_acc;
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

    # does this feature already exist?
    my $coords_q = "SELECT sfm.feature_id FROM seq_feat_mappings sfm, sequence_features f"
	. " WHERE sfm.seq_id = $seq_id"
	. " AND sfm.strand = \"$strand\""
	. " AND (($min >= sfm.feat_min AND $min <= sfm.feat_max)"
	. " OR ($max >= sfm.feat_min AND $max <= sfm.feat_max)"
	. " OR ($min <= sfm.feat_min AND $max >= sfm.feat_max))"
	. " AND f.feature_id=sfm.feature_id and f.feat_type = \"tRNA\"";
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
	my $id = $sacc . "-tRNA-$index";
	my $name = "tRNA-$aa-$AA{$aa}";

	my $featobj = new Bio::SeqFeature::Generic ( -start => $min,
						     -end => $max,
						     -strand => $strand,
						     -primary => 'tRNA',
						     -source_tag   => 'tRNAscan-SE',
						     -display_name => $id,
						     -score  => $score,
						     -tag => { product => $name, anticodon => $acodon });
	my $seqobj;
	my $fid = load_SeqFeature($dbh, $seq_id, $featobj, $seqobj, "0000253");
	
	my $annObj = Bio::Annotation::Collection->new();
	$annObj->add_Annotation("product",
				Bio::Annotation::SimpleValue->new(-value=>$name));
	$annObj->add_Annotation("anticodon",
				Bio::Annotation::SimpleValue->new(-value=>$acodon));
	load_feature_annotations($dbh, $fid, $annObj, "tRNAscan-SE", 3);
    }
}
