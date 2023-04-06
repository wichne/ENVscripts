#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
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
my $line = <$in>; # skip header line
while ($line = <$in>) {
    next if ($line =~ /^\#/);
    chomp $line;
    my ($acc,
	$fasta,
	$scaffold,
	$gene_pos,
	$start,
	$end,
	$strand,
	$rank,
	$kegg_id,
	$kegg_desc, # can be used as-is
	$pept_id, # merops
	$pept_family,
	$pept_hit, # use this, but trim off the '[pept_id] - ' at the beginning
	$pept_RBH, # true or false
	$pept_perid,
	$pept_bits,
	$pept_evalue,
	$pfam_hit,
	$CAZy_hit,
	$vogdb_hit,
	$vogdb_categories,
	$heme_reg_count) = split("\t", $line);

    my $seq_id = &get_seq_id_by_seq_accession($dbh, $scaffold);
	if (! $seq_id) { die "couldn't find seq_id from '$scaffold'\n"}

	my $feature_id = &get_feature_id_by_coords($dbh, $seq_id, "CDS", 
											$strand, $start, $end);
    if (! $feature_id) {
		warn "Couldn't get feature-id from $scaffold/$start-$end($strand)\n";
		next;
    }
    if ($kegg_desc) {
		my $source = "KEGG";
		my $rank = 10;
		my @desc = split(";", $kegg_desc);
		foreach my $desc (@desc) {
			my @ec;
			if ($kegg_desc =~ /(.*)\s+\[EC:(.*)\]/) {
				$desc = $1;
				my $ec = $2;
				@ec = split(/\s/, $ec);
			}
			$insert->execute($feature_id, $desc, "KEGG", "product");
			foreach my $ec(@ec) {
				$insert->execute($feature_id, $ec, "KEGG", "EC_number");
			}
			print "\tKEGG $desc\t" . join("\t", @ec) . "\n";
		}
    }

    if ($pept_hit) {
		my @desc = split(";", $pept_hit);
		foreach my $desc (@desc) {
			$desc =~ s/^$pept_id \- //;
			$desc =~ s/\(.*//;
			$insert->execute($feature_id, $desc, "MEROPS", "product");
			print "\tMEROPS $desc\n";
		}
    }

    if ($vogdb_hit) {
		my @desc = split(";", $vogdb_hit);
		foreach my $desc (@desc) {
			$desc =~ s/^\S+\s//;
			$desc =~ s/;.*$//;
			$insert->execute($feature_id, $desc, "VOGDB", "product");
			print "\tVOGDB $desc\n";
		}
    }

    if ($pfam_hit) {
		my @desc=split(";", $pfam_hit);
		foreach my $desc (@desc) {
			print "\tPFAM $desc\n";
		}
    }

    if ($CAZy_hit) {
		my @desc=split(";", $CAZy_hit);
		foreach my $desc(@desc) {
			print "\tCAZy $desc\n";
		}
    }
}
