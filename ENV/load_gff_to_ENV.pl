#!/usr/bin/perl

use strict;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Bio::DB::Fasta;
use Getopt::Std;
use lib $ENV{ENVSCRIPTS};
#use ProkGene;
use ENV;

my $opts = {};
&getopts('D:u:p:g:f:s:P:i:l:UN:h', $opts);
my $help = "load_gff_to_ENV.pl -D db -u user -p password -g gff_file -f fasta_file -s source -P acc_prefix -i set_id -N set_name -l locus_prefix -U \n"
    . " -U update coords and annotation for all features and insert new features (default is to only load new features)\n"
    . " -P prefix for accession (ie, GI, GP, etc)\n"
    . " -l locus prefix\n"
    . "\nSequences must be loaded before this program is run.\n";

if($opts->{h}) { print $help; exit() }
# if -U is not set, program will only load new information
# if -U is set, it will update coords and annotation for existing proteins and insert new proteins

my ($dbh, $gff_file, $fasta_file, $source,
    $locus_prefix, $prefix, $set_id, $update) = &handle_options($opts);
#        get the max existing locus tag
my ($max_locus_tag);
my $BYACC = {};
my $BYLOC = {};
if ($set_id) {
    my $max_locus_tag = max_locus_tag($dbh, $set_id);
    my $by_acc_q = "SELECT fx.feature_id, fx.accession from feature_accessions fx, seq_feat_mappings sfm, seq_set_link l WHERE set_id=$set_id and sfm.seq_id=l.seq_id and fx.feature_id=sfm.feature_id";
    my $by_acc_r = $dbh->selectall_arrayref($by_acc_q);
    foreach my $row (@$by_acc_r) {
	$BYACC->{$row->[1]} = $row->[0];
    }
    my $by_loc_q = "SELECT feature_id, feat_min, feat_max, l.seq_id, strand FROM seq_feat_mappings sfm, seq_set_link l WHERE set_id=$set_id and sfm.seq_id=l.seq_id";
    my $by_loc_r = $dbh->selectall_arrayref($by_loc_q);
    foreach my $row (@$by_loc_r) {
	$BYLOC->{$row->[3]}->{$row->[4]}->{$row->[1]} = $row->[0];
	$BYLOC->{$row->[3]}->{$row->[4]}->{$row->[2]} = $row->[0];
    }
}
#my $fadb = Bio::DB::Fasta->new($fasta_file, (-reindex => 1));
#my $FT = get_feat_types($dbh);

my $locus_index = 1;

my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);

my $SEQ = {};

 FEATURE:while (my $feature = $gffio->next_feature()) {
     # Look to see if there is a sequence accession (seq_id).
     # If there is, see if it's in the database.
     # If so, set the seq_id for the feature insert.
     # If not, skip insertion (with warning)
     my ($seqobj, $seq_id);
     my $seq_acc = $feature->seq_id;
     if ($seq_acc) {
	 if (defined $SEQ->{$seq_acc}) { $seqobj = $SEQ->{$seq_acc} }
	 else {
	     my @acc = split/\|/, $seq_acc;
	     #my @acc = split /\_/, $seq_acc;
	     foreach my $a (@acc) {
		 $seqobj = get_sequence_by_accession($dbh, $a, $set_id);
		 if ($seqobj) {
		     $SEQ->{$seq_acc} = $seqobj;
		     last;
		 }
	     }
	     if (! $seqobj) {
		 $seqobj = get_sequence_by_accession($dbh, $seq_acc, $set_id);
		 if (! defined $seqobj) {
		     print STDERR "Sequence '$seq_acc' is not in database. Skipping features associated with this accession\n";
		     next FEATURE;
		 } else {
		     $SEQ->{$seq_acc} = $seqobj;
		 }
	     }
	 }
	 
	 $seq_id = $seqobj->id;
	 
	 # set_id is global for the program, so this query should only be run once
	 if (! $set_id) {
	     $set_id = seq_id_to_set_id($dbh, $seq_id);
	     $max_locus_tag = max_locus_tag($dbh, $set_id);
	     my $by_acc_q = "SELECT fx.feature_id, fx.accession from feature_accessions fx, seq_feat_mappings sfm, seq_set_link l WHERE set_id=$set_id and sfm.seq_id=l.seq_id and fx.feature_id=sfm.feature_id";
	     my $by_acc_r = $dbh->selectall_arrayref($by_acc_q);
	     foreach my $row (@$by_acc_r) {
		 $BYACC->{$row->[1]} = $row->[0];
	     }
	     my $by_loc_q = "SELECT feature_id, feat_min, feat_max, l.seq_id, strand FROM seq_feat_mappings sfm, seq_set_link l WHERE set_id=$set_id and sfm.seq_id=l.seq_id";
	     my $by_loc_r = $dbh->selectall_arrayref($by_loc_q);
	     foreach my $row (@$by_loc_r) {
		 $BYLOC->{$row->[3]}->{$row->[4]}->{$row->[1]} = $row->[0];
		 $BYLOC->{$row->[3]}->{$row->[4]}->{$row->[2]} = $row->[0];
	     }
	 }

	 my ($max_locus_prefix, $max_locus_index);
	 if ($max_locus_tag) {
	     print STDERR "Setting max_locus_tag to $max_locus_tag\n";
	     if ($max_locus_tag =~ /^(.*\D)(\d+)(\.\d+)?$/) {
		 $max_locus_prefix = $1;
		 $max_locus_index = $2;
	     } else { warn "Can't parse max_locus_tag: '$max_locus_tag'\n"; }
	     if ($locus_prefix ne $max_locus_prefix) { warn "Locus prefix provided '$locus_prefix' is different from existing max locus_tag '$max_locus_tag'." }
	     else { $locus_index = $max_locus_index }
	 }
	 

     } else {
	 print STDERR "No sequence accession in row. Skipping...\n";
	 next FEATURE;
     }

     # grab data from object
     my $data = {};
     my ($min, $max, $strand) = ($feature->start, $feature->end, $feature->strand);
     if ($strand == "" or ! defined $strand) {
	 ($min, $max, $strand) = $min < $max ? ($min, $max, 1) : ($max, $min, -1);
     }
     my $feat_type = $feature->primary_tag();
     my @tags = $feature->get_all_tags();

     # Figure out what type of feature we have and handle it
     # BTW, the feat_type has to be in the INSDC feature table to be loaded
     # otherwise the process dies.

     if ($feat_type eq "source" || $feat_type eq "region") {
	 foreach my $tag (@tags) {
	     my @values = $feature->get_tag_values($tag);
	     # this currently isn't checking for loading redundant information
	     if ($tag eq "Is_circular" ||
		 $tag eq "Name" ||
		 $tag eq "genome" ||
		 $tag eq "mol_type") {
		 load_sequence_annotations($dbh, $seq_id, {$tag => $values[0]});
	     }
	 }
	 next FEATURE;
     }
     
     elsif ($feat_type eq "gene" ||
	 $feat_type eq "exon") { # we don't handle those right now
	 #warn "Skipping $feat_type ". $feature->display_name . "\n";
	 next FEATURE;
     } elsif ($feat_type =~ /repeat/) {
	 $feat_type = "repeat_region" if ($feat_type eq "repeat");
     } elsif ($feat_type eq "CDS" ||
	      $feat_type eq "ncRNA" ||
	      $feat_type eq "tRNA" ||
	      $feat_type eq "rRNA" ||
	      $feat_type eq "tmRNA"
	      ) {

	 my ($locus_tag, $max_partial, $min_partial);
	 foreach my $tag (@tags) {
	     my @values = $feature->get_tag_values($tag);
	     if ($tag eq "locus_tag") { $locus_tag = $values[0] }
	     
	     if ($tag eq "partial") {
		 $min_partial = $values[0] & 10 ? 1 : 0;
		 $max_partial = $values[0] & 01 ? 1 : 0;
	     }
	 }

	 my $fid;
	 # look for existing feature
	 # first by accession...
	 if ($locus_tag) {
	     if ($BYACC->{$locus_tag}) {
		 $fid = $BYACC->{$locus_tag};
	     }
	 }
	 # ...then by coords
	 if (! $fid) {
	     if (defined $BYLOC->{$seq_id}) {
		 if (defined $BYLOC->{$seq_id}->{$strand}) {
		     if (defined $BYLOC->{$seq_id}->{$strand}->{$max}) {
			 $fid = $BYLOC->{$seq_id}->{$strand}->{$max};
		     } elsif (defined $BYLOC->{$seq_id}->{$strand}->{$min}) {
			 $fid = $BYLOC->{$seq_id}->{$strand}->{$min};
		     }
		 }
	     }
	 }
	 
	 # if the feature is in the db, update the coords and annotation
	 if ($fid) {
	     if ($update) {
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
	     if (!$locus_tag) {
		 if ($feature->has_tag("ID")) {
		     my ($val) = $feature->get_tag_values("ID");
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
		     } elsif ($val =~ /^\d+$/) {
			 $locus_tag = sprintf "%s_%04d", ($locus_prefix, $val);
		     } else {
			 $locus_tag = $val;
		     }
		 } elsif ($locus_prefix && $locus_index) {
		     $locus_tag = sprintf "${locus_prefix}_%05d", $locus_index++;
		 }
		 $feature->set_attributes(-tag => {'locus_tag' => $locus_tag });
	     }
	     
	     # insert the information
	     warn "inserting new row for $locus_tag on $seq_id\n";
	     my $SO_term;
	     my $feat_id = load_SeqFeature($dbh, $seq_id, $feature, $seqobj, $SO_term, $source, $prefix);
	 }
     } else {
	 warn "Oooh! $feat_type is a new type. I need to learn how to handle it! ($gff_file)\n";
	 next FEATURE;
     }
}

sub handle_options {
    my $opts = shift;
    my $err_msg;
    
    my $dbh = &connect($opts);
    my $gff_file = $opts->{'g'} || {$err_msg .= "Need to specify gff file to load with -g\n"};
    my $fasta_file = $opts->{'f'};
    my $set_id = $opts->{'i'};
    #    if (! $fasta_file && ! $seq_id) {$err_msg .= "Need to specify fasta file of underlying metagenome sequence with -f or seq_id with -i\n"}
    if (! $set_id && $opts->{'N'}) {
	$set_id = set_name_to_id($dbh, $opts->{'N'});
	if (!$set_id) { $err_msg .= "Need to either provide a set_id with -i or a set_name  with -N\n";}
    }
    my $source = $opts->{'s'} || {$err_msg .= "Need to specify information source (Genbank, IMG, etc.) with -s\n"};
    my $prefix = $opts->{'P'} ? $opts->{'P'} : "";
    my $locus_prefix = $opts->{'l'};
    my $update = $opts->{'U'};
    
    if ($err_msg) { die $err_msg }

    return ($dbh, $gff_file, $fasta_file, $source,
	    $locus_prefix, $prefix, $set_id, $update);
}

sub get_feat_types {
    my $dbh = shift;
    my $q = 'select * from INSDC.feature_key';
    my $r = $dbh->selectall_hashref($q, 'feature_key');
    return $r;
}
