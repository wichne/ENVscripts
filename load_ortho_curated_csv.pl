#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;

my %arg;
&getopts('f:D:p:u:c:h', \%arg);

if ($arg{'h'}) {
    print "USAGE load_MR_xlsx.pl -f csvfile -D database -p dbpassword [ -u dbuser -c curator -h ]

-f path to the csv file that contains the curated annotation. There must be columns with the headers 'Product', and/or 'Gene', and/or 'EC' for this program to load anything. Also, there must be an accession or identifier by which the program can look up feature_id in the database.

-D the database to update
-u database user (optional) Default is current username.
-p the mysql password
-c if the user specified is uploading someone else's work, you can specify the curation author with this option

-h This help message
";
    exit();
}
my $user = $arg{'u'} ? $arg{'u'} : $ENV{USER};

my $dbh = connect(\%arg);

$| = 1;

open (my $in, $arg{f}) or die "Can't open input file $arg{f}: $!\n";
# the first row is headers, or should be.
# Look for a locus_tag field or accession field
# "Product"
# "Gene"
# "EC"
# some of these may not exist
my $first_line = <$in>;
if ($first_line !~ /Product/) { die "The first line needs to be a header line as specified in the (-h) help\n"; }
chomp $first_line;
my @cell = split/\t/, $first_line;
my ($product_col,
    $gene_col,
    $ec_col,
    @acc_col);
for (my $col=0; $col<@cell; $col++) {
    if ($cell[$col] =~ /^\s*product\s*/i) {
	$product_col = $col;
    } elsif ($cell[$col] =~ /^\s*gene\s*$/i) {
	$gene_col = $col;
    } elsif ($cell[$col] =~ /^\s*ec\s*/i) {
	$ec_col = $col;
    }
    else {
	push @acc_col, $col;
    }
}
	
# Now parse through the data rows, pulling the cols of interest
# and just printing them out for now.
my $lookup; # find an accession to lookup the feature_id
while (my $line = <$in>) {
    chomp $line;
    my @cell = split/\t/, $line;
    
    my ($description,
	$gene,
	$ec);
    my $product = $cell[$product_col];
    if ($product) {
	$description = $product;
    }
    my $gene = $cell[$gene_col];
    my $ec = $cell[$ec_col];
    
    ##  if there's no annotation values
    ## don't load anything:
    if (! ($description || $gene || $ec)) { next }
	    
    foreach my $col (@acc_col) {
	if ($cell[$col]) {
	    my $acc_value = $cell[$col];
	    my @accs = split/[\s\,]+/, $acc_value;
	    
	    foreach my $lookup (@accs) {
		my @x = split/\|/, $lookup;
		my $fid;
		foreach my $lu (@x) {
		    if ($lu =~ /C7867-\d{3}/) {next}
		    $fid = get_feature_id_by_accession($dbh, $lu);
		    last if ($fid);
		}
		if (! $fid) {
		    print STDERR "Couldn't find feature_id for $lookup\n";
		    next;
		}
		
		my $featref = get_features_by_feature_id($dbh, $fid);
		my $annotationref = $featref->{$fid}->{'annotation'}; # rank, source, qualifier
		# if there is existing curated annotation, set is_current=0

		my $AnnColl = Bio::Annotation::Collection->new();
		if (defined $annotationref->{2}) {
		    ## If there's an existing curated row with the same information,
		    ## preserve it, and don't load a duplicate row.
		    ## Otherwise set the is_current for the old annotation to 0
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
		    my $source = $arg{'c'} ? $arg{'c'} : $user;
		    load_feature_annotations($dbh, $fid, $AnnColl, $source, 2);
		} else {
		    $AnnColl->add_Annotation('product', Bio::Annotation::SimpleValue->new(-value => $description));
		    $AnnColl->add_Annotation('gene', Bio::Annotation::SimpleValue->new(-value => $gene));
		    $AnnColl->add_Annotation('EC_number', Bio::Annotation::SimpleValue->new(-value => $ec));
		    # load the current data
		    my $source = $arg{'c'} ? $arg{'c'} : $user;
		    load_feature_annotations($dbh, $fid, $AnnColl, $source, 2);
		    # print STDERR "\n";
		}
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
