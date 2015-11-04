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
while (my $line = <$IN>) {
    next if ($line =~ /^#/);
    chomp $line;
    my ($seq_acc,
	$source,
	$method,
	$start,
	$end,
	$score,
	$dir,
	$frame,
	$attributes) = split /\s+/, $line;

    my @accs = split/\|/, $seq_acc;
    my $sid;
    my $sacc;
    foreach my $ax(@accs) {
	$sid = get_seq_id_by_accession($dbh, $ax);
	$sacc = $ax;
	if ($sid) { last }
    }
    if (!$sid) { die " Couldn't get a sequence_id from the accession '$seq_acc'.\n" }

    my $strand = $dir eq "-" ? -1 : 1;

    my @atts = split(/\;/, $attributes);
    my ($id, $model);
    foreach my $at (@atts) {
	my ($key, $val) = split/\=/, $at;
	if ($key eq "rfam-id") {
	    $RFAM{$val}++;
	    my $index = $RFAM{$val};
	    $id = $sacc . "-" . $val . "-" . $index;
	} elsif ($key eq "rfam-acc") {
	    $model = $val;
	}
    }

    my ($type, $name, $gene, $so_term) = &annotation_from_Rfam($dbh, $model);

    my $featobj = new Bio::SeqFeature::Generic ( -start => $start,
						 -end => $end,
						 -strand => $strand,
						 -primary => $type,
						 -source_tag   => 'Rfam',
						 -display_name => $id,
						 -score  => $score );
    
    my $seqobj;
    my $fid = load_SeqFeature($dbh, $sid, $featobj, $seqobj, $so_term);

    my $annObj = Bio::Annotation::Collection->new();
    $annObj->add_Annotation("product",
			    Bio::Annotation::SimpleValue->new(-value=>$name));
    $annObj->add_Annotation("gene",
			    Bio::Annotation::SimpleValue->new(-value=>$gene)) if ($gene);

    delete_feature_annotations($dbh, $fid);
    load_feature_annotations($dbh, $fid, $annObj, "Rfam", 3);
}

sub annotation_from_Rfam {
    my $dbh = shift;
    my $model = shift;
    my $query = "SELECT feature_type, product, gene, SO_term from rfam.rfam where accession=\"$model\"";
    my $row = $dbh->selectrow_arrayref($query);
    return (@$row);
}
