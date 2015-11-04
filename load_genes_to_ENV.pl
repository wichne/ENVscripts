#!/usr/bin/perl

use lib "/home/nels329/devel";
use ENV;
use strict;
use Getopt::Std;
use Bio::SeqFeature;

my $args = {};
&getopts('i:d:u:p:', $args);

my $input_file=$args->{i};

open (my $IN, $input_file) or die "Can't open $input_file: $!\n";

while (my $line = <$IN>) {
    chomp $line;
    my ($id,
	$end5,
	$end3,
	$dir,
	$annotation) = split /\t/, $line, 5;
    my $sid = get_seq_id_by_accession($dbh, $id);
    if (!$sid) { die " Couldn't get a sequence_id from the accession '$id'.\n" }

    my $featobj = ;
    my $fid = load_SeqFeature($dbh, $sid, $featobj, $seqobj, $SO_term);

    my $annobj = Bio::Annotation::Collection->new();
    foreach my $key("product", "gene", "EC_number") {
	if (defined $ref->{$key}) {
	    my $sv = Bio::Annotation::SimpleValue->new(-value=>$ref->{$key});
	    $annObj->add_Annotation($key, $sv);
	}
    }
    delete_feature_annotations($dbh, $feature_id);
    load_feature_annotations($dbh, $feature_id, $annObj);
}
