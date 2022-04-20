#!/usr/bin/perl

use strict;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Bio::DB::Fasta;
use Getopt::Std;
use URI::Escape;
use lib $ENV{ENVSCRIPTS};
#use ProkGene;
use ENV;

my $opts = {};
&getopts('D:u:p:g:f:s:P:i:l:', $opts);

my ($dbh, $gff_file, $fasta_file, $source, $locus_prefix, $prefix, $set_id) = &handle_options($opts);
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
	    $seqobj = get_sequence_by_accession($dbh, $a, $set_id);
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
    } elsif ($feat_type eq "mRNA" || $feat_type eq "exon" || $feat_type eq "gene") { next FEATURE }

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
	    # make a PNNL id, which is the contig acc plus the CDS index
	    if ($seq_acc =~ /(.*\_scaf)\_(\d+)$/) {
		my $scaf_prefix = $1;
		my $seqidx = $2;
		if ($val =~ /^\d+\_(\d+)$/) {
		    $source = "PNNL";
		    my $cdsidx = $1;
		    $locus_tag = sprintf "%s_%d_%d", ($scaf_prefix, $seqidx, $cdsidx);
		    print "Making a PNNL id: $locus_tag\n";
		    $feature->set_attributes(-tag=>{'locus_tag' => $locus_tag });
		}
	    } else {
		$locus_tag = $val;
		$feature->set_attributes(-tag => {'locus_tag' => $locus_tag });
	    }
	}
    }
    print STDERR "Doing $feat_type $locus_tag...\n";
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

    
    # if the feature is in the db, update the coords and annotation and add the accession (if necessary)
    if (@$feat_r) {
	print STDERR "\tFound matching feature in db: " . $feat_r->[0] . "\n";
	delete_feature_annotations($dbh, $feat_r->[0], {'source' => $source, 'ann_rank' => 10});
	my $feat_ann_i = "INSERT feature_annotations"
	    . " (feature_id, data_type_id, value, source, ann_rank)"
	    . " VALUES (?, ?, ?, ?, ?)";
	my $ann_h = $dbh->prepare($feat_ann_i);

	if (grep /^product$/, @tags) {
	    my $product = uri_unescape([$feature->get_tag_values('product')]->[0]);
	    $ann_h->execute($feat_r->[0], 66, $product, $source, 10);
	} elsif (grep /^Name$/, @tags) {
	    my $product = uri_unescape([$feature->get_tag_values('Name')]->[0]);
            $ann_h->execute($feat_r->[0], 66, $product, $source, 10);
        }

	if (grep /^EC$/, @tags) {
	    foreach my $ec (@{$feature->get_tag_values('EC')}) {
		$ann_h->execute($feat_r->[0], 1, $ec, $source, 10);
	    }
	}

	if (grep /^sso$/i, @tags) { # KBase SEED sequence ortholog
	    $ann_h->execute($feat_r->[0], 23, [$feature->get_tag_values('sso')]->[0], $source, 10);
	}

	# insert the accession
	my $feat_acc_i = "INSERT feature_accessions"
	    . " (feature_id, accession, source, prefix)"
	    . " VALUES (?, ?, ?, ?)";
	my $acc_h = $dbh->prepare($feat_acc_i);
	
	my @accs = split(/\|/, $locus_tag);
	for my $acc (@accs) {
	    if ($acc ~~ ["gb", "RefSeq", "PNNL", "RAST", "gp", "gi", "fig"]) { next }
	    else { $acc_h->execute($feat_r->[0], $acc, $source, $prefix) }
	}
	

       # update seq_feat_link

    } else {
	# insert the information
	my $SO_term;
	my $feat_id = load_SeqFeature($dbh, $seq_id, $feature, $seqobj, $SO_term, $source, $prefix);
 	my $feat_ann_i = "INSERT feature_annotations"
	    . " (feature_id, data_type_id, value, source, ann_rank)"
	    . " VALUES (?, ?, ?, ?, ?)";
	my $ann_h = $dbh->prepare($feat_ann_i);

	if (grep /^product$/, @tags) {
	    my $product = uri_unescape([$feature->get_tag_values('product')]->[0]);
	    $ann_h->execute($feat_r->[0], 66, $product, $source, 10);
	} elsif (grep /^Name$/, @tags) {
	    my $product = uri_unescape([$feature->get_tag_values('Name')]->[0]);
            $ann_h->execute($feat_r->[0], 66, $product, $source, 10);
        }

	if (grep /^EC$/, @tags) {
	    foreach my $ec (@{$feature->get_tag_values('EC')}) {
		$ann_h->execute($feat_r->[0], 1, $ec, $source, 10);
	    }
	}

	if (grep /^sso$/i, @tags) { # KBase SEED sequence ortholog
	    $ann_h->execute($feat_r->[0], 23, [$feature->get_tag_values('sso')]->[0], $source, 10);
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
    my $set_id = $opts->{'i'};

    if ($err_msg) { die $err_msg }

    return ($dbh, $gff_file, $fasta_file, $source, $locus_prefix, $prefix, $set_id);
}

sub get_feat_types {
    my $dbh = shift;
    my $q = 'select * from INSDC.feature_key';
    my $r = $dbh->selectall_hashref($q, 'feature_key');
    return $r;
}
