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
&getopts('D:u:p:g:s:i:N:hk:', $opts);
my $help = "load_gff_to_ENV.pl -D db -u user -p password -g gff_file -f fasta_file -s source -P acc_prefix -i set_id -N set_name -l locus_prefix -U \n"
    . " -U update coords and annotation for all features and insert new features (default is to only load new features)\n"
    . " -P prefix for accession (ie, GI, GP, etc)\n"
    . " -l locus prefix\n"
    . "\nSequences must be loaded before this program is run.\n";

if($opts->{h}) { print $help; exit() }

my ($dbh, $gff_file, $source, $acc_key,
    $set_id) = &handle_options($opts);

# store the existing features in the database by acc and by location for faster retrieval later
my $BYACC = {};
my $BYLOC = {};
if ($set_id) {
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

my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);

my $SEQ = {}

 FEATURE:while (my $feature = $gffio->next_feature()) {
     # Look to see if there is a sequence accession (seq_id).
     # If there is, see if it's in the database.
     # If so, set the seq_id for the feature insert.
     # If not, skip insertion (with warning)
     my ($seqobj, $seq_id);
     my $seq_acc = $feature->seq_id;
     if ($seq_acc) {
	 if ($SEQ{$seq_acc}) { $seq_id = $SEQ{$seq_acc} }
	 else {
	     my @acc = split/\|/, $seq_acc;
	     foreach my $a (@acc) {
		 $seq_id = get_seq_id_by_seq_accession($dbh, $a);
		 if ($seq_id) {
		     $SEQ{$seq_acc} = $seq_id;
		     last;
		 }
	     }
	 }
	 if (! defined $seq_id) {
	     print STDERR "Sequence '$seq_acc' is not in database. Skipping features associated with this accession\n";
	     next FEATURE;
	 }
	 
	 # set_id is global for the program, so this query should only be run once
	 if (! $set_id) {
	     $set_id = seq_id_to_set_id($dbh, $seq_id);
	     
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

     } else {
	 print STDERR "No sequence accession in row? Skipping...\n";
	 next FEATURE;
     }

     # grab data from object
     my $data = {};
     my $feat_type = $feature->primary_tag();

     # Figure out what type of feature we have and handle it
     # BTW, the feat_type has to be in the INSDC feature table to be loaded
     # otherwise the process dies.

     if ($feat_type eq "source" || $feat_type eq "region") {
	 next FEATURE;
     }
     
     elsif ($feat_type eq "gene" ||
	 $feat_type eq "exon") { # we don't handle those right now
	 #warn "Skipping $feat_type ". $feature->display_name . "\n";
	 next FEATURE;
     } elsif ($feat_type =~ /repeat/) {
	 next FEATURE;
     } elsif ($feat_type eq "CDS" ||
	      $feat_type eq "ncRNA" ||
	      $feat_type eq "tRNA" ||
	      $feat_type eq "rRNA" ||
	      $feat_type eq "tmRNA"
	      ) {

	 my ($locus_tag, $new_acc);

	 if ($feature->has_tag($acc_key)) {
	     my @values = $feature->get_tag_values($acc_key);
	     if (@values) { $new_acc = $values[0]; }
	 }
	 if ($feature->has_tag("locus_tag")) {
	     my @values = $feature->get_tag_values("locus_tag");
	     if (@values) { $locus_tag = $values[0] }
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
	     my ($min, $max, $strand) = ($feature->start, $feature->end, $feature->strand);
	     if ($strand == "" or ! defined $strand) {
		 ($min, $max, $strand) = $min < $max ? ($min, $max, 1) : ($max, $min, -1);
	     }
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
	 
	 # if the feature is in the db, insert the accessions
	 if ($fid && $new_acc) {
	     insert_feature_accessions($dbh, $fid, $new_acc, $source);
	 } else {
	     print "There was a problem: '$seq_id' : '$locus_tag' -> '$fid' -> '$new_acc'\n";
	 }
     } else {
	 next FEATURE;
     }
}

sub handle_options {
    my $opts = shift;
    my $err_msg;
    
    my $dbh = &connect($opts);
    my $gff_file = $opts->{'g'} || {$err_msg .= "Need to specify gff file to load with -g\n"};
    my $set_id = $opts->{'i'};
    #    if (! $fasta_file && ! $seq_id) {$err_msg .= "Need to specify fasta file of underlying metagenome sequence with -f or seq_id with -i\n"}
    if (! $set_id && $opts->{'N'}) {
	$set_id = set_name_to_id($dbh, $opts->{'N'});
	if (!$set_id) { $err_msg .= "Need to either provide a set_id with -i or a set_name  with -N\n";}
    }
    my $source = $opts->{'s'} || {$err_msg .= "Need to specify information source (Genbank, IMG, etc.) with -s\n"};
    my $acc_key = $opts->{'k'} || die "Need to provide a metadata key to load as an accession";
	
    if ($err_msg) { die $err_msg }

    return ($dbh, $gff_file, $source, $acc_key,
	    $set_id);
}

sub get_feat_types {
    my $dbh = shift;
    my $q = 'select * from INSDC.feature_key';
    my $r = $dbh->selectall_hashref($q, 'feature_key');
    return $r;
}
