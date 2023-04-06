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
&getopts('D:u:p:a:f:s:P:', $opts);
# if -U is not set, program will only load new information
# if -U is set, it will update coords and annotation for existing proteins and insert new proteins

my ($dbh, $annotation_file, $source) = &handle_options($opts);

#query  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description     Preferred_name  GOs     EC      KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE   KEGG_TC CAZy    BiGG_Reaction   PFAMs
#Contig11747308_6        1379270.AUXF01000006_gene80     4.67e-55        184.0   COG0847@1|root,COG0847@2|Bacteria,1ZTNH@142182|Gemmatimonadetes 142182|Gemmatimonadetes L       EXOIII  -       -       2.7.7.7 ko:K02342       ko00230,ko00240,ko01100,ko03030,ko03430,ko03440,map00230,map00240,map01100,map03030,map03430,map03440   M00260  R00375,R00376,R00377,R00378     RC02795 ko00000,ko00001,ko00002,ko01000,ko03032,ko03400 -       -       -       RNase_T

open(my $in, $annotation_file) or die "Can't open $annotation_file as input: $!\n";
while (my $line = <$in>) {
    next if ($line =~ /^#/);
    chomp $line;
    my ($query,
	$seed_ortholog,
	$evalue,
	$score,
	$eggNOG_OGs,
	$max_annot_lvl,
	$COG_category,
	$description,
	$preferred_name,
	$GOs,
	$EC,
	$KEGG_ko,
	$KEGG_pathway,
	$KEGG_module,
	$KEGG_reaction,
	$KEGG_rclass,
	$BRITE,
	$KEGG_TC,
	$CAZy,
	$BiGG_reaction,
	$PFAMs) = split(/\t/, $line);

    # look for existing feature
    # first by accession...
    my $feat_r = [];
    if ($query) {
	my $acc_q = "SELECT feature_id from feature_accessions where accession = \"$query\"";
	$feat_r = $dbh->selectcol_arrayref($acc_q);
    }
     
    # if the feature is in the db, update the coords and annotation
    if (@$feat_r) {
	my $feat_id = $feat_r->[0];
	if (@$feat_r > 1) {
	    warn "Feature $query has multiple feat_ids associated: " . join(", ", @$feat_r) . ". Using the first.\n";
	}
	# what are we going to load? eggNOG_OGs, description (or preferred name?), GOs, ECs, KEGG_ko, TC, CAZy,
	delete_feature_annotations($dbh, $feat_id, {'source' => 'emapper'});
	delete_feature_evidence($dbh, $feat_id, {'program' => 'emapper'});
	my $feature_evidence_i = "INSERT feature_evidence"
	    . " (feature_id, feat_min, feat_max, ev_type, ev_accession, ev_date,"
	    . " program, score)"
	    . " VALUES($feat_id, 1, 1, ?, ?, NOW(), ?, 1)";
	my $feature_annotation_i = "INSERT feature_annotations"
	    . " (feature_id, data_type_id, value, edit_by, date_edit, ann_rank, source)"
	    . " SELECT $feat_id, d.id, ?, USER(), NOW(), ?, ?"
	    . " FROM INSDC.qualifier d"
	    . " WHERE d.qualifier=?";
	
	# eggNOG
	my @NOGs = split(",", $eggNOG_OGs);
	my %SEEN;
	foreach my $nog(@NOGs) {
	    my @p = split('@', $nog);
	    if ($p[0] =~ /^COG\d{4}/) {
		if (!$SEEN{$p[0]}) {
		    my $prod = $dbh->selectcol_arrayref("SELECT description from egad.COG where accession=\"$p[0]\"");
		    $description=$prod->[0];
		    $dbh->do($feature_evidence_i, {},
			     ('eggNOG', $p[0], 'emapper')) or warn $DBI::errstr;
		    $SEEN{$p[0]}=1;
		}
	    }
	}

	# description
	if ($description ne "-") {
	    $dbh->do($feature_annotation_i, {},
		     ($description, 0, 'emapper', 'product')) or warn $DBI::errstr;
	}
	# preferred_name = gene_sym
	if ($preferred_name ne "-") {
	    $dbh->do($feature_annotation_i, {},
		     ($preferred_name, 0, 'emapper', 'gene')) or warn $DBI::errstr;
	}
	# GO terms
	if ($GOs ne "-") {
	    my @gos = split(",", $GOs);
	    foreach my $go(@gos) {
		$dbh->do($feature_annotation_i, {},
			 ($go, 0, 'emapper', 'db_xref')) or warn $DBI::errstr;
	    }
	}
	# EC terms
	if ($EC ne "-") {
	    my @ecs = split(",", $EC);
	    foreach my $ec(@ecs) {
		$dbh->do($feature_annotation_i, {},
			 ($ec, 0, 'emapper', 'EC_number')) or warn $DBI::errstr;
	    }
	}
	# KO
	if ($KEGG_ko ne "-") {
	    my @kos = split(",", $KEGG_ko);
	    foreach my $ko(@kos) {
		$dbh->do($feature_evidence_i, {},
			 ('KO', $ko, 'emapper')) or warn $DBI::errstr;
	    }
	}
	# TC terms
	if ($KEGG_TC ne "-") {
	    my @tcs = split(",", $KEGG_TC);
	    foreach my $tc(@tcs) {
		$dbh->do($feature_annotation_i, {},
			 ("TC|".$tc, 0, 'emapper', 'db_xref')) or warn $DBI::errstr;
	    }
	}
	# CAZy
	if ($CAZy ne "-") {
	    my @cazys = split(",", $CAZy);
	    foreach my $c(@cazys) {
		$dbh->do($feature_annotation_i, {},
			 ("CAZy|".$c, 0, 'emapper', 'db_xref')) or warn $DBI::errstr;
	    }
	}
    }
}

sub handle_options {
    my $opts = shift;
    my $err_msg;
    
    my $dbh = &connect($opts);
    my $input_file = $opts->{'a'} || {$err_msg .= "Need to specify input file to load with -a\n"};
    my $source = $opts->{'s'};
    
    if ($err_msg) { die $err_msg }

    return ($dbh, $input_file, $source);
}

sub get_feat_types {
    my $dbh = shift;
    my $q = 'select * from INSDC.feature_key';
    my $r = $dbh->selectall_hashref($q, 'feature_key');
    return $r;
}
