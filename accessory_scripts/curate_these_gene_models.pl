#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use Spreadsheet::ParseXLSX;
use Spreadsheet::ParseExcel;

my %arg;
&getopts('f:D:p:u:c:i:h', \%arg);

if ($arg{'h'}) {
    print "USAGE curate_these_gene_models.pl -f xlsxfile -D database -p dbpassword -i set_id [ -u dbuser -c curator -h ]

-f path to the xlsx file that contains the curated annotation. There must be columns with the headers 'Curated product descriptor', and/or 'Gene', and/or 'Curated EC' for this program to load anything. Also, there must be an accession or identifier by which the program can look up feature_id in the database.

-D the database to update
-u database user (optional) Default is current username.
-p the mysql password
-i set_id
-c if the user specified is uploading someone else's work, you can specify the curation author with this option

-h This help message
";
    exit();
}
my $user = $arg{'u'} ? $arg{'u'} : $ENV{USER};

my $dbh = connect(\%arg);

$| = 1;
my $parser = Spreadsheet::ParseXLSX->new;
my $wkbk = $parser->parse($arg{f});
my $set_id = $arg{i};

for my $wksht ($wkbk->worksheets() ) {
    # Could parse the 'drop' table, but instead
    # will keep track of genes in db vs genes in 'pro'
    # table and set f.is_current=0 for those that are not
    # in 'pro' table.
    print STDERR "Found worksheet ", $wksht->get_name, "\n";
    if ($wksht->get_name =~ /pro$/) {
	my %CURRENT;
	my ( $row_min, $row_max ) = $wksht->row_range();
	my ( $col_min, $col_max ) = $wksht->col_range();
	
	# the first row is headers, or should be.
	# Look for a locus_tag field or accession field
	# "Curated Product Descriptor"
	# "Gene"
	# "Curated EC"
	# some of these may not exist
	my ($product_col,
	    $gene_col,
	    $ec_col,
	    $img_col,
	    $rast_col,
	    $locus_col);
	for my $col ($col_min .. $col_max) { 
	    my $cell = $wksht->get_cell($row_min, $col);
	    if (! defined $cell) { print STDERR "No cell: $row_min, $col\n"; next; }
	    if ($cell->value =~ /curated product descriptor/i) {
		$product_col = $col;
	    } elsif ($cell->value =~ /^\s*gene\s*$/i) {
		$gene_col = $col;
	    } elsif ($cell->value =~ /curated ec/i) {
		$ec_col = $col;
	    } elsif ($cell->value =~ /rast acc/i) {
		$rast_col = $col;
	    } elsif ($cell->value =~ /img acc/i) {
		$img_col = $col;
	    } elsif ($cell->value =~ /locus[\s\_]tag/i) {
		$locus_col = $col;
	    }
	}
	
	# Now parse through the data rows, pulling the cols of interest
	# and just printing them out for now.
	for my $row ( $row_min+1 .. $row_max ) {
	    my $lookup; # find an accession to lookup the feature_id
	    if (defined $wksht->get_cell($row, $locus_col) &&
		$wksht->get_cell($row, $locus_col)->value) {
		$lookup = $wksht->get_cell($row, $locus_col)->value;
	    } elsif (defined $wksht->get_cell($row, $img_col) &&
		     $wksht->get_cell($row, $img_col)->value) {
		$lookup = $wksht->get_cell($row, $img_col)->value
	    } elsif (defined $wksht->get_cell($row, $rast_col) &&
		     $wksht->get_cell($row, $rast_col)->value) {
		$lookup = $wksht->get_cell($row, $rast_col)->value; }
	    # print STDERR "Found $lookup...";
	    my $fid = get_feature_id_by_accession($dbh, $lookup);
	    if (! $fid) {
		print STDERR "Couldn't find feature_id for $lookup\n";
		next;
	    }
	    $CURRENT{$fid} = 1;
	}
	# Here is where I want to compare what's in the database 
	# against what's in CURRENT and is_current=0 those that
	# aren't there
	my $all_fid_q = "SELECT f.feature_id from sequence_features f, seq_feat_mappings m, seq_set_link l"
	    . " WHERE l.set_id=$set_id"
	    . " AND m.seq_id=l.seq_id"
	    . " AND f.feature_id=m.feature_id"
	    . " AND f.is_current=1"
	    . " AND f.feat_type=\"CDS\"";
	my $fids = $dbh->selectcol_arrayref($all_fid_q);
	my $is_c_q = "UPDATE sequence_features set is_current=0 where feature_id=?";
	my $sth = $dbh->prepare($is_c_q);
	foreach my $f (@$fids) {
	    if (! defined $CURRENT{$f}) {
		print STDERR "$f\t?\n";
# Let's not do this now. It looks like more work needs to be done on this issue.
#		$sth->execute($f);
	    }
	}
    } elsif ($wksht->get_name =~ /drop$/i) {
	my ( $row_min, $row_max ) = $wksht->row_range();
	my ( $col_min, $col_max ) = $wksht->col_range();
	my @locus_cols;
	for my $col ($col_min .. $col_max) { 
	    my $cell = $wksht->get_cell($row_min, $col);
	    if (! defined $cell) {
#		print STDERR "No cell: $row_min, $col\n";
		next;
	    }
	    if ($cell->value =~ /locus[\s\_]tag/i) {
#		print STDERR "FOUND locus : $col\n";
		push @locus_cols, $col;
	    }
	}
	for my $row ( $row_min+1 .. $row_max ) {
	    my $lookup; # find an accession to lookup the feature_id
	    foreach my $col (@locus_cols) {
#		print STDERR "Looking in $row, $col\n";
		if (defined $wksht->get_cell($row, $col) &&
		    $wksht->get_cell($row, $col)->value) {
		    $lookup = $wksht->get_cell($row, $col)->value;
		    last;
		}
	    }
	    print "$lookup\tDEL\n";
	}
    }
}
