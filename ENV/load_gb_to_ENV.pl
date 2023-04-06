#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use vars qw ($dbh %FEAT_NAME);
$| = 1;
my $args = {};
&getopts('D:u:p:i:F:P:s:S:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : 'genbank';
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $locus_prefix = $args->{'P'};
my $set_id = $args->{'s'}; # or die "Need to provide set_id with -s\n";
my $source = $args->{'S'};

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
    my $seq_q = "SELECT s.seq_id, sequence, seq_accession"
	. " FROM sequences s, sequence_accessions a, seq_set_link l"
	. " WHERE a.seq_accession = \"" . $seqrec->display_name
	. "\" AND l.seq_id=a.seq_id"
	. " AND l.set_id=$set_id"
	. " AND s.seq_id=a.seq_id";
    my $seq_q_r = $dbh->selectall_hashref($seq_q, 'seq_id');
    if (scalar(keys %$seq_q_r) > 1) { die "ERROR:: Sequence accession " . $seqrec->display_name . " has multiple seq_ids in database for set_id $set_id. Dying.\n";}
    foreach my $sid(keys %$seq_q_r) {
	$seq_id = $sid;
	$dbseqobj = Bio::Seq->new(-seq => $seq_q_r->{$sid}->{'sequence'},
				  -id  => $sid,
				  -accession => $seq_q_r->{$sid}->{'seq_accession'});
    }

    # Is the first Genbank record feature the sequence?
    if (! $dbseqobj) {
	warn "Couldn't find a seq_id for " . $seqrec->display_name . " so I'm inserting it.\n";
        # insert the sequence into assembly, stan and asmbl_data
	$seq_id = &load_sequence_SeqFeature($dbh, $seqrec);
	print STDERR "Sequence " . $seqrec->display_name . " loaded: $seq_id\n";

	# if there's no set_id, load a set
	if (! $set_id) {
	    warn "Loading a set: $filename\n";
	    my $name = $filename;
	    $name =~ s/\.[^\.]+$//;
	    $set_id = load_sequence_sets($dbh, $name, $seqrec->desc);
	    if (! $set_id) { die "Can't do a sequence_sets insertion\n"; }
	    else { warn "\t$set_id loaded.\n"; }
	}
	link_seq_to_set($dbh, $set_id, $seq_id);
#	print STDERR "\tand linked to $set_id\n";
	$dbseqobj = get_sequence_by_seq_id($dbh, $seq_id);
    }

    # Otherwise skip this entry with a warning
#    else {
 #       warn "Sequence " . $seqrec->accession . " is not in database or in the file(?). Skipping loading features associated with this accession\n";
#	next;
 #   }

    #============== The Features ============================#
    foreach my $feature (@features) {
	my $feat_type = $feature->primary_tag();
	# Don't load the sequence as a feature
	if ($feat_type eq "source") { next }
	# load all but 'gene' features (except gene features that are pseudogenes - load those).
	if ($feat_type eq "gene" && !($feature->has_tag('pseudo'))) { next }
	#if ($feat_type ne "CDS") { print "Not dealing with $feat_type\n"; next; }
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

	my $feat_r;
	if (! $locus_tag) { $locus_tag = "NEW" . $feature->primary_tag . "$min"; }
	else {
	    # look for existing feature
	    # first by accession...
	    my $acc_q = "SELECT feature_id from feature_accessions where accession = \"$locus_tag\"";
	    $feat_r = $dbh->selectcol_arrayref($acc_q);
	}
	
	# ...then by coords
	if (! defined $feat_r || ! @$feat_r) {
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
		if ($feature->annotation->get_all_annotation_keys) {
		    &delete_feature_annotations($dbh, $fid, {'source' => $source});
		    &load_feature_annotations($dbh, $fid, $feature->annotation, $source, 10);
		}

	    }
	}
	# Otherwise must insert the new information
	else {
	    my $feat_id = &load_SeqFeature($dbh, $seq_id, $feature, $dbseqobj, undef, $source);
	    print "inserted $feat_id " . $locus_tag . " at " . $feature->start . " / " . $feature->end . "\n";
	    if ($feature->annotation->get_all_annotation_keys) {
		&load_feature_annotations($dbh, $feat_id, $feature->annotation, $source, 10);
	    }
	    $DONE{$feat_id} = 1;
	}
    }
}

print "\n";
exit();

