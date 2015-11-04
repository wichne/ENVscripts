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
    print "USAGE load_MR_xlsx.pl -f xlsxfile -D database -p dbpassword -i set_id [ -u dbuser -c curator -h ]

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
my $set_id = $arg{i} or die "Need to provide set_id with -i\n";

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
	    }
	    
	    if ($description && $ftcolor eq "#000000") {
		# print STDERR "updating.\n";
		my $featref = get_features_by_feature_id($dbh, $fid);
		my $annotationref = $featref->{$fid}->{'annotation'}; # rank, source, qualifier
		# if there is existing curated annotation, set is_current=0
		if (defined $annotationref->{2}) {
		    while (my ($user, $qualref) = each %{$annotationref->{2}}) {
			foreach my $val (@{$qualref->{'product'}}) {
			    if ($val ne $description) {
				&set_annotation_not_current($dbh, $fid, 2, $user, 66, $val);
				$AnnColl->add_Annotation('product', Bio::Annotation::SimpleValue->new(-value => $description));
			    }
				
			}
			foreach my $val (@{$qualref->{'gene'}}) {
			    if ($val ne $gene) {
				&set_annotation_not_current($dbh, $fid, 2, $user, 35, $val);
				$AnnColl->add_Annotation('gene', Bio::Annotation::SimpleValue->new(-value => $gene));
			    }
			}
			foreach my $val (@{$qualref->{'EC_number'}}) {
			    if ($val ne $ec) {
				&set_annotation_not_current($dbh, $fid, 2, $user, 1, $val);
				$AnnColl->add_Annotation('EC_number', Bio::Annotation::SimpleValue->new(-value => $ec));
			    }
			}
		    }
		} else {
		    $AnnColl->add_Annotation('product', Bio::Annotation::SimpleValue->new(-value => $description));
		    $AnnColl->add_Annotation('gene', Bio::Annotation::SimpleValue->new(-value => $gene));
		    $AnnColl->add_Annotation('EC_number', Bio::Annotation::SimpleValue->new(-value => $ec));
		}
		# load the current data
		my $source = $arg{'c'} ? $arg{'c'} : $user;
		load_feature_annotations($dbh, $fid, $AnnColl, $source, 2);
	    } else {
		#print STDERR "Not updating.\n";
	    }
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
		print STDERR "$f is not listed in spreadsheet.\n";
# Let's not do this now. It looks like more work needs to be done on this issue.
#		$sth->execute($f);
	    }
	}
    }
}

sub set_annotation_not_current {
    my ($dbh,
	$fid,
	$rank,
	$source,
	$data_type_id,
	$val) = @_;
    my $q = "UPDATE feature_annotations SET is_current=0"
	. " WHERE feature_id=$fid"
	. " AND rank = $rank"
	. " AND source = \"$source\""
	. " AND data_type_id = $data_type_id"
	. " AND value = \"$val\""; 
    $dbh->do($q); 
}
