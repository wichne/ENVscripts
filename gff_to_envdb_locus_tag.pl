#!/usr/bin/perl

#Modified original gff_to_envdb.pl Jmo 08/25/2015 to preserve prodigal gene calls as accession of laminar metagenome, see line 74

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
	    if ($val =~ /^(\d+)\_(\d+)$/) {
		# This is a prodigal ID. Change it to desired format
		my $featidx = $2;
		my $scaf= $feature->seq_id; #grab seq_id which contains scaffold name
		$locus_tag= sprintf "%s_%06d", ($scaf, $featidx);
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
    if (@$feat_r) {
	print "Accession $locus_tag is already in db. skipping...\n";
	next;
    }

    # ...then by coords
    else {
	my $coords_q = "SELECT feature_id FROM seq_feat_mappings"
	    . " WHERE seq_id = $seq_id"
	    . " AND strand = \"$strand\""
	    . " AND (feat_min = $min OR feat_max = $max)";
	$feat_r = $dbh->selectcol_arrayref($coords_q);

	if (!@$feat_r) {
	    print "Couldn't find a feature using coords $min/$max. Skipping...\n";
	} else {
	    insert_feature_accessions($dbh, $feat_r->[0], $acc, $source, $prefix);
	}
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
