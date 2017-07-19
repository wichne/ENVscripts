#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use ENV;
use vars qw ($dbh %FEAT_NAME);
$| = 1;
my $args = {};
&getopts('D:u:p:i:F:P:s:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $locus_prefix = $args->{'P'};
my $set_id = $args->{'s'};

my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
$dbh = DBI->connect("dbi:mysql:host=$host;database=$db",
		       $user, $pswd, {'AutoCommit'=> 0});
# Keep track of which feature_ids we update in case of coordinate conflicts/overlapping genes
my %DONE;

while (my $seqrec = $in->next_seq) {
    my @features = $seqrec->get_SeqFeatures();

    # ============= The sequence =========================#
    # if the sequence is already in the database, retrieve the seq_id, else load it
    my ($dbseqobj, $seq_id);
    $dbseqobj = get_sequence_by_accession($dbh, $seqrec->accession);
    if ($dbseqobj) { $seq_id = $dbseqobj->id }

    # Is the first Genbank record feature the sequence?
    elsif ($features[0]->primary_tag eq "source") {
        # insert the sequence into assembly, stan and asmbl_data
	my $seqObj = shift @features;
	$seq_id = &load_sequence_SeqFeature($dbh, $seqObj);
	print STDERR "Sequence " . $seqObj->display_id . " loaded: $seq_id\n";
	link_seq_to_set($dbh, $set_id, $seq_id);
#	print STDERR "\tand linked to $set_id\n";
	$dbseqobj = get_sequence_by_seq_id($dbh, $seq_id);
    }

    # Otherwise skip this entry with a warning
    else {
        warn "Sequence " . $seqrec->accession . " is not in database or in the file(?). Skipping loading features associated with this accession\n";
	next;
    }

    #============== The Features ============================#
    foreach my $feature (@features) {
	my $feat_type = $feature->primary_tag();
	# Don't load the sequence as a feature
	if ($feat_type eq "source") { next }
	# load all but 'gene' features (except gene features that are pseudogenes - load those).
#	if ($feat_type eq "gene" && !($feature->has_tag('pseudo'))) { next }
	if ($feat_type ne "CDS") { print "Not dealing with $feat_type\n"; next; }
	# Get the location of the feature
	my $data = {};
	my ($min, $max, $strand) = ($feature->start, $feature->end, $feature->strand);
	if ($strand == "" or ! defined $strand) {
	    ($min, $max, $strand) = $min < $max ? ($min, $max, 1) : ($max, $min, -1);
	}
	my $FT_coords = $feature->location->to_FTstring;

	# Get an accession from the tags
	my @tags = $feature->get_all_tags();
	my ($locus_tag, $max_partial, $min_partial);
	foreach my $tag (@tags) {
	    my @values = $feature->get_tag_values($tag);
	    if ($tag eq "locus_tag") { $locus_tag = $values[0] }
	    if ($tag eq "partial") {
		$min_partial = $values[0] & 10 ? 1 : 0;
		$max_partial = $values[0] & 01 ? 1 : 0;
	    }
	}
	if (! $locus_tag) { $locus_tag = "NEW" . $feature->primary_tag . "$min"; }

	# look for existing feature
	# first by accession...
	my $acc_q = "SELECT feature_id from feature_accessions where accession = \"$locus_tag\"";
	my $feat_r = $dbh->selectcol_arrayref($acc_q);
	# ...then by coords
	if (! @$feat_r) {
	    my $coords_q = "SELECT sfm.feature_id FROM seq_feat_mappings sfm, sequence_features f"
		. " WHERE seq_id = $seq_id"
		. " AND strand = \"$strand\""
#		. " AND (($min >= sfm.feat_min AND $min <= sfm.feat_max)"
#		. " OR ($max >= sfm.feat_min AND $max <= sfm.feat_max)"
#		. " OR ($min <= sfm.feat_min AND $max >= sfm.feat_max))"
		. " AND ($min = sfm.feat_min OR $max = sfm.feat_max)"
		. " AND f.feature_id=sfm.feature_id and f.feat_type = \"$feat_type\"";
	    $feat_r = $dbh->selectcol_arrayref($coords_q);
	}
	
	# if the feature is in the db
	if (@$feat_r) {
	    foreach my $fid (@$feat_r) {
		if ($DONE{$fid}) { print "feature $fid has already been updated. Must have conflicting coordinates. Please check this one manually.\n"; next; }
		my $featr = get_features_by_feature_id($dbh, $fid);
		# Only update if the information differs
		$DONE{$fid} = 1;
		if ($featr->{$fid}->{'location'}->{$seq_id}->{'feat_min'} != $min ||
		    $featr->{$fid}->{'location'}->{$seq_id}->{'feat_max'} != $max ||
		    $featr->{$fid}->{'location'}->{$seq_id}->{'strand'} != $strand ||
		    $featr->{$fid}->{'location'}->{$seq_id}->{'min_partial'} != $min_partial ||
		    $featr->{$fid}->{'location'}->{$seq_id}->{'max_partial'} != $max_partial) {
		    print "Update $fid $locus_tag from " . 
			$featr->{$fid}->{'location'}->{$seq_id}->{'feat_min'} . " / " .
			$featr->{$fid}->{'location'}->{$seq_id}->{'feat_max'} .
			" to $min/$max ($FT_coords)\n";
		    update_feature_mapping($dbh, $seq_id, $fid, {'feat_min'    => $min,
								 'feat_max'    => $max,
								 'strand'      => $strand,
		                                                 'translation_coords' => $FT_coords,
								 'min_partial' => $min_partial,
								 'max_partial' => $max_partial});
		}
	    }
	}
	# Otherwise must insert the new information
	else {
	    my $feat_id = &load_SeqFeature($dbh, $seq_id, $feature, $dbseqobj, undef, "JGI");
#	    my $feat_id;
	    print "inserted $feat_id " . $locus_tag . " at " . $feature->start . " / " . $feature->end . "\n";
	    $DONE{$feat_id} = 1;
	}
    }
}

print "\n";
exit();

