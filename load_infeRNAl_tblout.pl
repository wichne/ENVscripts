#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Bio::SeqFeature::Generic;
use Getopt::Std;


my %arg;
&getopts('D:u:p:i:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
if (! $infile) { die "No input file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}


my $feature_annotation_i = "INSERT feature_annotations"
    . " (feature_id, data_type_id, value, edit_by, date_edit, ann_rank, source)"
    . " SELECT ?, d.id, ?, USER(), NOW(), 10, ?"
    . " FROM INSDC.qualifier d"
    . " WHERE d.qualifier=?";
my $insert = $dbh->prepare($feature_annotation_i);

open my $in, $infile or die "Can't open $infile: $!\n";
my %IDX;
while (my $line = <$in>) {
    next if ($line =~ /^\#/);
    chomp $line;
    my ($idx,
	$type,
	$accession,
	$scaffold,
	$dash,
	$clan,
	$model,
	$mdl_lo,
	$mdl_hi,
	$start, # this is not low
	$end,
	$strand,
	$trunc,
	$pass,
	$gc,
	$bias,
	$score,
	$evalue,
	$inc,
	$olp,
	$anyidx,
    $afrct1,
    $afrct2,
    $winidx,
    $wfrct1,
    $wfrct2,
	$description) = split(/\s+/, $line, 27);

    if ($description eq "-") { $description = $type }

    if ($type eq "5S_rRNA") { $type = "rRNA" }
    elsif ($type eq "SSU_rRNA_bacteria") { $type = "rRNA" }
    elsif ($type eq "LSU_rRNA_bacteria") { $type = "rRNA" }
    elsif ($type eq "tRNA") { $type = "tRNA" }
    elsif ($type eq "tmRNA") { $type = "tmRNA" }
    else { $type = "misc_RNA" }

    $IDX{$type}++;

    if ($strand eq "+") { $strand = 1;} elsif ($strand eq "-") { $strand = -1 }
    my $seq_id = &get_seq_id_by_seq_accession($dbh, $scaffold);
	if (! $seq_id) { die "couldn't find seq_id from '$scaffold'\n"}

	my $feature_id = &get_feature_id_by_coords($dbh, $seq_id, $type, 
											$strand, $start, $end);
    if (! $feature_id) {
        my ($min, $max) = sort {$a<=>$b} ($start, $end);
        my $feature_acc = sprintf "%s-%s-%05i", ($scaffold,$type,$IDX{$type});
        my $featObj = new Bio::SeqFeature::Generic ( -start => $min,
                                -end => $max,
                                -strand => $strand,
                                -primary => $type,
                                -source_tag   => 'infeRNAl',
                                -display_name => $feature_acc,
                                -score  => $score );
        $featObj->add_tag_value("product",$description);
        
        my $seqobj;
        my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, undef, undef, "infeRNAl");
        if (! $feat_id) { die "Failed to load feature" }
    }
}
