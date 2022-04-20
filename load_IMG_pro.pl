#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;
use DBI;
use Spreadsheet::ParseXLSX;
use Spreadsheet::ParseExcel;

my %arg;
&getopts('i:D:p:h', \%arg);
my $source = "IMG";
my $rank = 9;

if ($arg{'h'}) {
    print STDERR "load_IMG_pro.pl -i [inputfile] -D [ENVdatabase] -p [dbpswd]

This program loads annotation from the excel-format IMG .pro.xls file
into the feature_annotations table in an ENV-schema database.

Datatypes loaded: Product_name, gene, EC

The program DOES delete or existing data in the database (for any gene loaded).

";
    exit;
}

my $db = $arg{'D'};
my $pswd = $arg{'p'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $ENV{USER}, $pswd);

my $parser = Spreadsheet::ParseXLSX->new;
my $wkbk = $parser->parse($arg{i});
for my $wksht ($wkbk->worksheets() ) {
    my ( $row_min, $row_max ) = $wksht->row_range();
    my ( $col_min, $col_max ) = $wksht->col_range();
    my ($product_col,
	$gene_col,
	$ec_col,
	$acc_col);
    for my $col ($col_min .. $col_max) { 
	my $cell = $wksht->get_cell($row_min, $col);
	if (! defined $cell) { print STDERR "No cell: $row_min, $col\n"; next; }
	if ($cell->value =~ /Product Name/i) {
	    $product_col = $col;
	} elsif ($cell->value =~ /Gene Symbol/i) {
	    $gene_col = $col;
	} elsif ($cell->value =~ /Enzymes/i) {
	    $ec_col = $col;
	} elsif ($cell->value =~ /Locus Tag/i) {
	    $acc_col = $col;
	}
    }
    for my $row ( $row_min+1 .. $row_max ) {
	my $lookup; # find an accession to lookup the feature_id
	if (defined $wksht->get_cell($row, $acc_col) &&
	    $wksht->get_cell($row, $acc_col)->value) {
	    $lookup = $wksht->get_cell($row, $acc_col)->value;
	}
	my $fid = get_feature_id_by_accession($dbh, $lookup);
	if (! $fid) {
	    print STDERR "Couldn't find feature_id for $lookup\n";
	    next;
	}
	my $AnnColl = Bio::Annotation::Collection->new();
	my ($description,
	    $ftcolor,
	    $gene,
	    $ec);
	my $product = $wksht->get_cell($row, $product_col);
	if (defined $product) {
	    $description = $product->value;
	    # if the product text is red, don't load it.
	    my $format = $product->get_format;
	    my $font = $format->{Font};
	    $ftcolor = $font->{Color};
	}
	
	if ($gene_col && defined $wksht->get_cell($row, $gene_col)) {
	    $gene = $wksht->get_cell($row, $gene_col)->value ;
	}
	if ($ec_col && defined $wksht->get_cell($row, $ec_col)) {
	    $ec = $wksht->get_cell($row, $ec_col)->value;
	    $ec =~ s/EC\://g;
	}
	$AnnColl->add_Annotation('product', Bio::Annotation::SimpleValue->new(-value => $description));
	$AnnColl->add_Annotation('gene', Bio::Annotation::SimpleValue->new(-value => $gene));
	my @ec = split(/\s+/, $ec);
	foreach my $e(@ec) {
	    $AnnColl->add_Annotation('EC_number', Bio::Annotation::SimpleValue->new(-value => $e));
	}

	# Delete any old annotation
	delete_feature_annotations($dbh, $fid, {'source' => 'IMG'});
	# load the current data
	load_feature_annotations($dbh, $fid, $AnnColl, $source, $rank);
    }
}
