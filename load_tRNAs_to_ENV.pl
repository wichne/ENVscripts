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
while (my $line = <$IN>) {
    chomp $line;
    my ($seq_acc,
	$index,
	$end5,
	$end3,
	$aa,
	$acodon,
	$intron_end5,
	$intron_end3,
	$score) = split /\s+/, $line;
    if ($index !~ /\d+/) { next }
    my @accs = split/\|/, $seq_acc;
    my $sid;
    my $sacc;
    foreach my $ax(@accs) {
	$sid = get_seq_id_by_accession($dbh, $ax);
	$sacc = $ax;
	if ($sid) { last }
    }
    if (!$sid) { die " Couldn't get a sequence_id from the accession '$seq_acc'.\n" }

    my ($start, $end, $strand) = $end5 < $end3 ? ($end5, $end3, 1) : ($end3, $end5, -1);
    $AA{$aa}++;
    my $id = $sacc . "-tRNA-$index";
    my $name = "tRNA-$aa-$AA{$aa}";

    my $featobj = new Bio::SeqFeature::Generic ( -start => $start,
						 -end => $end,
						 -strand => $strand,
						 -primary => 'tRNA',
						 -source_tag   => 'tRNAscan-SE',
						 -display_name => $id,
						 -score  => $score,
						 -tag => { product => $name, anticodon => $acodon });
    my $seqobj;
    my $fid = load_SeqFeature($dbh, $sid, $featobj, $seqobj, "0000253");

    my $annObj = Bio::Annotation::Collection->new();
    $annObj->add_Annotation("product",
			    Bio::Annotation::SimpleValue->new(-value=>$name));
    $annObj->add_Annotation("anticodon",
			    Bio::Annotation::SimpleValue->new(-value=>$acodon));

    delete_feature_annotations($dbh, $fid);
    load_feature_annotations($dbh, $fid, $annObj, "tRNAscan-SE", 3);
}
