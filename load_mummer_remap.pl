#!/usr/bin/perl
use lib $ENV{'SCRIPTS'};
use ENV;
use Getopt::Std;
use strict;

my %args;
&getopts('D:p:i:s:', \%args);
my $file = $args{'i'} || die "provide input with -i\n";
my $source = $args{'s'} || die "provide source with -s\n";
open (my $map, $file);

my $dbh = connect(\%args);

# 0 - fasta file
# 1 - seq_acc
# 2 - feat_type
# 3 - locus_tag
# 4 - low-coord
# 5 - hi-coord
# 6 - strand
# 7 - ref fasta file
# 8 - ref seq_acc
# 9 - ref feat_type
# 10 - ref locus_tag
# 11 - ref low coord
# 12 - ref high coord
# 13 - ref strand

while (my $line = <$map>) {
    chomp $line;
    my @f = split(/\t/, $line);
# rows with problem mappings start with "#"
    if ($line =~ /^#/) {
	if ($f[10] && $f[3] !~ /NEW/) {
	    my $ref_fid = &get_feature_id_by_accession($dbh, $f[10]);
	    my $qry_fid = &get_feature_id_by_accession($dbh, $f[3]);
	    if (!($ref_fid && $qry_fid)) { warn "Hey! $f[10] $f[3]\n"; next; }
# if the row starts with a "#", then there is a significant difference in the gene model. The assumption is that 
# the query genes will be more complete than the ref genes.
# OR if there is no ref feature info
# Insert a sequence_features row using info in $_[2] and calculate a protein sequence (if possible) from 
# $_[1], $_[4], $_[5], $_[6].
# Insert a seq_feat_mappings row using $_[1], $_[4], $_[5], $_[6].
# Insert a feature_accessions row using $_[3]. If this is a # row, insert a row for $_[10], too.
	    my ($prefix, $accession) = $f[10] =~ /\|/ ? split(/\|/, $f[10], 2) : ("", $f[10]);
	    my $sx_q = "insert feature_accessions (feature_id, source, accession)"
		. " values ($qry_fid, \"$source\", \"$accession\")";
#	    $dbh->do($sx_q) || print "$sx_q\n";
	}
    } elsif ($f[10] && $f[3] !~ /NEW/) {
# if a row doesn't start with "#", grab the feature_id from the accession in field $_[10],
# and insert a seq_feat_mappings row with that feature_id, and the location info in fields $_[1], $_[4], $_[5], $_[6]
# Then insert a feature_accessions row using the info in $_[3]. If $_[3] is "NEW", increment the max locus_tag.
	my $ref_fid = &get_feature_id_by_accession($dbh, $f[10]);
	my $qry_fid = &get_feature_id_by_accession($dbh, $f[3]);
	if (! ($ref_fid && $qry_fid)) { warn "Ho! $f[10] $ref_fid $f[3] $qry_fid\n"; next; }
	my $sfm_q = "update seq_feat_mappings set feature_id=$ref_fid where feature_id=$qry_fid";
	my $sx_q = "update feature_accessions set feature_id=$ref_fid where feature_id=$qry_fid and (source=\"locus_tag\" or source = \"RAST\")";
#	$dbh->do($sfm_q) || print "$sfm_q\n";
#	$dbh->do($sx_q) || print "$sx_q\n";
    }
}

exit();

