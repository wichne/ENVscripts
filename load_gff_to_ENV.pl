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
#my $fadb = Bio::DB::Fasta->new($fasta_file, (-reindex => 1));
#my $FT = get_feat_types($dbh);

my $locus_index = 1;

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
	 if (! $seq_id) { $seqobj = get_sequence_by_accession($dbh, $seq_acc, $set_id) }
	 if (! defined $seqobj) {
	     print STDERR "Sequence '$seq_acc' is not in database. Skipping features associated with this accession\n";
	     next FEATURE;
	 }
	 $seq_id = $seqobj->id;
	 if (! $set_id) {
	     $set_id = seq_id_to_set_id($dbh, $seq_id);
        }
#        get the max existing locus tag
	 my $max_locus_tag = max_locus_tag($dbh, $set_id);
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
     if ($feat_type eq "gene" ||
	 $feat_type eq "exon") { # we don't handle those right now
	 warn "Skipping $feat_type ". $feature->display_name . "\n";
	 next FEATURE;
     }
     $feat_type = "repeat_region" if ($feat_type eq "repeat");

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
     }

     # look for existing feature
     # first by accession...
     my $feat_r = [];
     if ($locus_tag) {
	 my $acc_q = "SELECT feature_id from feature_accessions where accession = \"$locus_tag\"";
	 $feat_r = $dbh->selectcol_arrayref($acc_q);
     }
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
	 if ($update) {
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
	 }
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
