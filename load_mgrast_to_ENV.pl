#!/usr/bin/perl

use DBI;
use strict;
use Getopt::Std;

my %arg;
&getopts('c:a:o:g:D:u:p:s:',\%arg);

my $db = $arg{'D'} or die "Need to provide -D database\n";
my $user = $arg{'u'} ? $arg{'u'} : $ENV{USER};
my $pswd = $arg{'p'} or die "Need to provide -p pswd\n";
my $SOURCE = $arg{'s'} ? $arg{'s'} : 7; # SEED by default

my $dbh = DBI->connect("dbi:mysql:server=localhost;db=$db", $user, $pswd);

my %DATA;

open my $ann, $arg{'a'} or die "Provide valid -a annotation_file (.superblat.expand.protein): $!\n";
while (my $line = <$ann>) {
    chomp $line;
    my ($M5_id,
	$q_id,
	$per_id,
	$length,
	$evalue,
	$function_id,
	$organism_id,
	$source) = split(/\t/, $line);

    if ($source != $SOURCE) { next }
    $DATA{$q_id}->{'function_id'} = $function_id;
}
close $ann;

if ($arg{'o'}) {
    open my $ont, $arg{'o'} or die "Invalid -o ontology_file (.superblat.expand.ontology): $!\n";
    while (my $line = <$ont>) {
	chomp; $line;
	my ($M5_id,
	    $q_id,
	    $per_id,
	    $length,
	    $evalue,
	    $function_id,
	    $ontology_id,
	    $source) = split(/\t/, $line);
	push @{$DATA{$q_id}->{'ontology'}->{$source}}, $ontology_id;
    }
    close $ont;
}
my %CLUSTER;
open my $clust, $arg{'c'} or die "Provide valid -c cluster_file (.cluster.aa90.mapping): $!\n";
while (my $line = <$clust>) {
    chomp $line;
    my ($clust_id,
	$key_id,
	$others,
	$per_id) = split(/\t/, $line);
    my @members = ($key_id, split(/\,/, $others));
    foreach my $m (@members) {
	$CLUSTER{$m} = $clust_id;
    }
}
close $clust;
open my $gen, $arg{'g'} or die "Provide valid -g gene_file (.genecalling.coding.faa headers): $!\n";
while (my $line = <$gen>) {
    chomp $line;
    if ($line =~ /^>?((.*)\_\[cov\=([\d\.]+)\]\_(\d+)\_(\d+)\_([\+\-]))/) {
	my ($racc,
	    $seq_acc,
	    $cov,
	    $min,
	    $max,
	    $dir) = ($1, $2, $3, $4, $5, $6);
	my $strand = $dir eq "+" ? 1 : $dir eq "-" ? -1 : 0;

	# Insert the gene object
	my $feat_i = "INSERT INTO sequence_features"
	    . " (feat_type_id, SO_term, inserted_by, date_inserted)"
	    . " VALUES (6, 316, USER(), NOW())";
	$dbh->do($feat_i);
	my $feat_id = $dbh->last_insert_id(undef, undef, undef, undef);
	if (!$feat_id) { die "what happened to $racc feat_id?\n"; }

	# Insert the location on the contig
	my $mapping_i = "INSERT INTO seq_feat_mappings (seq_id, feature_id, feat_min, feat_max, strand)"
	    . " SELECT seq_id, $feat_id, $min, $max, \"$strand\" FROM sequence_accessions"
	    . " WHERE seq_accession = \"$seq_acc\"";
	$dbh->do($mapping_i) || warn "Problem with query: $mapping_i";

	my $ann_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, edit_by, date_edit, rank) VALUES ($feat_id,?,?, USER(), NOW(), ?)";

	#Insert the MGRAST accession
	$dbh->do($ann_i, {}, (23, "MGRAST:$racc", undef));

	#Insert the aa90 cluster id, if any
	$dbh->do($ann_i, {}, (23, "aa90:$CLUSTER{$racc}", undef)) if ($CLUSTER{$racc});

	#Get the functional annotation
	my $product = "protein of unknown function";
	my $func_id = $DATA{$racc}->{'function_id'} ? $DATA{$racc}->{'function_id'} : 
	    $DATA{$CLUSTER{$racc}}->{'function_id'};
	if ($func_id) {
	    my $func_q = "SELECT function from M5NR.function where id=$func_id";
	    ($product) = $dbh->selectrow_array($func_q);
	}
	if (! $product) { warn "Problem with $racc func_id $func_id: no function found." }

	#Find and insert EC numbers
	while ($product =~ /(\d+(\.(\d+|\-){3})/g) {
	    $dbh->do($ann_i, {}, (1, $1, undef));
	}

	#Insert the product descriptor
	$dbh->do($ann_i, {}, (66, $product, undef));

	#Insert ontology terms
	foreach my $ont (keys %{$DATA{$racc}->{'ontology'}}) {
	    foreach my $term (@{$DATA{$racc}->{'ontology'}->{$ont}}) {
		my $q = "INSERT INTO feature_annotations (feature_id, data_type_id, value, edit_by, date_edit)"
		    . " SELECT $feat_id, 23, concat(source,':',accession), USER(), NOW() FROM M5NR.ontology o, M5NR.source s"
		    . " WHERE o.id=$term AND s.id = o.source_id";
		$dbh->do($q);
	    }
	}
    }
}
