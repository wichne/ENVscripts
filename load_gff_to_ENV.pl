#!/usr/bin/perl

use strict;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Bio::DB::Fasta;
use Getopt::Std;
use lib $ENV{SCRIPTS};
#use ProkGene;
use ENV;

my $opts = {};
&getopts('D:u:p:g:f:s:P:i:l:', $opts);

my ($dbh, $gff_file, $fasta_file, $source, $locus_prefix, $prefix) = &handle_options($opts);
#my $fadb = Bio::DB::Fasta->new($fasta_file, (-reindex => 1));
#my $FT = get_feat_types($dbh);

my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);

FEATURE:while (my $feature = $gffio->next_feature()) {
    # Look to see if there is a sequence accession (seq_id).
    # If there is, see if it's in the database.
    # If so, set the seq_id for the feature insert.
    # If not, skip insertion (with warning)
    my ($seqobj,$seq_id);
    my $seq_acc = $feature->seq_id;
    if ($seq_acc) {
	my @acc = split/\|/, $seq_acc;
	#my @acc = split /\_/, $seq_acc;
	foreach my $a (@acc) {
	    $seqobj = get_sequence_by_accession($dbh, $a);
	    last if $seqobj;
	}
	$seq_id = $seqobj->id;
	if (! defined $seqobj) {
	    print STDERR "Sequence '$seq_acc' is not in database. Skipping features associated with this accession\n";
	    next FEATURE;
	}
    } else {
	print STDERR "No sequence accession in row. Skipping...\n";
	next FEATURE;
    }

    my $feat_type = $feature->primary_tag();
    if ($feat_type eq "source") {
#	my $asmbl_id_q = "SELECT a.asmbl_id FROM assembly a, stan s"
#	    . " WHERE length(sequence) = " . $feature->end()
#	    . " AND iscurrent = 1 AND s.asmbl_id=a.asmbl_id";
#	my $res = $sgc->dbh->selectall_arrayref($asmbl_id_q);
#	if (@$res > 1) { die "Too many current asmbl_ids with length " . $feature->end . "\n" }
#	elsif (@$res == 0) { die "No asmbl_ids of length " . $feature->end() . "\n" }
#	else { $asmbl_id = $res->[0]->[0] }
	next FEATURE;
    }

    # grab data from object
    my $data = {};
    my ($min, $max, $strand) = ($feature->start, $feature->end, $feature->strand);
    if ($strand == "" or ! defined $strand) {
	($min, $max, $strand) = $min < $max ? ($min, $max, 1) : ($max, $min, -1);
    }
    my @tags = $feature->get_all_tags();
    my ($locus_tag, $max_partial, $min_partial);
    foreach my $tag (@tags) {
	my @values = $feature->get_tag_values($tag);
	if ($tag eq "locus_tag") { $locus_tag = $values[0] }
	if ($tag eq "partial") {
	    $min_partial = $values[0] & 10 ? 1 : 0;
	    $max_partial = $values[0] & 01 ? 1 : 0;
	}
	if ($tag eq "ID" && !$locus_tag) {
	    my ($val) = $feature->get_tag_values($tag);
	    if (! $locus_prefix) { 
		print "No locus tag for feature ID=$val. Please enter a prefix so a meaningful locus_tag can be created: ";
		my $entry = <STDIN>;
		chomp $entry;
		$locus_prefix = $entry;
	    }
	    if ($val =~ /^(\d+)\_(\d+)$/) {
		# This is a prodigal ID. Change it to another format
#		$feature->seq_id =~ /(\d+)/;
		my $seqidx = $1;
		$val =~ /\_(\d+)/;
		my $featidx = $1;
		$locus_tag = sprintf "%s_%05d_%04d", ($locus_prefix, $seqidx, $featidx);
		$feature->set_attributes(-tag => {'locus_tag' => $locus_tag });
	    } else {
		$locus_tag = $val;
		$feature->set_attributes(-tag => {'locus_tag' => $locus_tag });
	    }
	}
    }

    # look for existing feature
    # first by accession...
	
    my $acc_q = "SELECT feature_id from feature_accessions where accession = \"$locus_tag\"";
    my $feat_r = $dbh->selectcol_arrayref($acc_q);

    # ...then by coords
    if (! @$feat_r) {
	my $coords_q = "SELECT feature_id FROM seq_feat_mappings"
	    . " WHERE seq_id = $seq_id"
	    . " AND strand = \"$strand\""
	    . " AND (feat_min = $min OR feat_max = $max)";
	$feat_r = $dbh->selectcol_arrayref($coords_q);
    }

    
    # if the feature is in the db, update the coords and annotation
    if (@$feat_r) {
	foreach my $fid (@$feat_r) {
	    my $featr = get_features_by_feature_id($dbh, $fid);
	    if ($featr->{$fid}->{'location'}->{$seq_id}->{'feat_min'} != $min ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'feat_max'} != $max ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'strand'} != $strand ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'min_partial'} != $min_partial ||
		$featr->{$fid}->{'location'}->{$seq_id}->{'max_partial'} != $max_partial) {
		update_feature_mapping($dbh, $seq_id, $fid, {'feat_min'    => $min,
							     'feat_max'    => $max,
							     'strand'      => $strand,
							     'min_partial' => $min_partial,
							     'max_partial' => $max_partial});
	    }
	}
	# update feature_annotations
    } else {
	# insert the information
	warn "inserting new row for $locus_tag on $seq_id\n";
	my $SO_term;
#	my $feat_id = load_SeqFeature($dbh, $seq_id, $feature, $seqobj, $SO_term, $source, $prefix);
    }
}

sub handle_options {
    my $opts = shift;
    my $err_msg;

    my $dbh = &connect($opts);
    my $gff_file = $opts->{'g'} || {$err_msg .= "Need to specify gff file to load with -g\n"};
    my $fasta_file = $opts->{'f'};
    my $seq_id = $opts->{'i'};
#    if (! $fasta_file && ! $seq_id) {$err_msg .= "Need to specify fasta file of underlying metagenome sequence with -f or seq_id with -i\n"}
    my $source = $opts->{'s'} || {$err_msg .= "Need to specify information source (Genbank, IMG, etc.) with -s\n"};
    my $prefix = $opts->{'P'} ? $opts->{'P'} : "";
    my $locus_prefix = $opts->{'l'};

    if ($err_msg) { die $err_msg }

    return ($dbh, $gff_file, $fasta_file, $source, $locus_prefix, $prefix);
}

sub get_feat_types {
    my $dbh = shift;
    my $q = 'select * from INSDC.feature_key';
    my $r = $dbh->selectall_hashref($q, 'feature_key');
    return $r;
}
