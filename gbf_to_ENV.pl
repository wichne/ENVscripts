#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:f:s:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $source = $args->{'s'} or die "Need to provide information source with -s\n";

# crack open the Genbank file:
my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

# crack open the peptide fasta file
my $pep = $args->{'f'};
my $dbo;
if ($pep) {
    $dbo = Bio::DB::Fasta->new($pep);
}

$dbh = DBI->connect("dbi:mysql:host=$host;database=$db",
		    $user, $pswd, {'AutoCommit'=> 0});

my $set_id; # = load_sequence_sets($dbh, $filename, "", 0);

while (my $seqo = $in->next_seq) {
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence into assembly, stan and asmbl_data

    
    my $seq_id = &get_seq_id_by_seq_accession($dbh, $seqo->display_name);
    if (! $seq_id) {
	$seq_id = &load_sequence_SeqFeature($dbh, $seqo);
	if (! $seq_id) { die "Failed loading sequence " . $seqo->display_name . "\n"; }
	else { print STDERR "Sequence " . $seqo->display_id . " loaded: $seq_id\n"; }
	link_seq_to_set($dbh, $set_id, $seq_id);
    }
    
    my ($start, $end, $locus_tag);
    foreach my $featObj (@features) {
	$start = $featObj->start;
	$end = $featObj->end;
	my $feat_id = &get_feature_id_by_coords($dbh, $seq_id, $start, $end);
	if ($featObj->primary_tag eq "gene") {
	    if ($featObj->has_tag('locus_tag')) {
		$locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
	    }
	    if ($featObj->has_tag('pseudo') && $pep) {
		my $aa = $dbo->get_Seq_by_id($featObj->display_name);
		if ($feat_id) {
		    print STDERR "MSG: $locus_tag already loaded\n";
		    # update sequence_features.product field
		    &update_product($dbh, $feat_id, $aa);
		} else {
		    $featObj->add_tag_value('translation', $aa);
		    $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqo, undef, "JGI");
		    print "inserted $feat_id " . $featObj->get_tag_values('locus_tag') . " at " . $featObj->start . " / " . $featObj->end . "\n";
		}
	    }
	    next;
	}
	
	# Only reach this part if the feature is not a gene (eg, CDS)
	# This sets the locus_tag for the CDS to the locus_tag harvested from the gene
	if ($featObj->start eq $start && $featObj->end eq $end && $locus_tag) {
	    $featObj->set_attributes(-tag => {'locus_tag' => $locus_tag });
	}
	
	# This is a hack for DRAM gbk files (sigh)
	if (! $locus_tag && $featObj->has_tag('gene')) {
	    $locus_tag = [$featObj->get_tag_values('gene')]->[0];
	    $featObj->set_attributes(-tag => {'locus_tag' => $locus_tag});
	}
	
	if ($feat_id) {
	    # the feature is already loaded, only update information
	    if ($locus_tag) {
		my $this_feat_id = &get_feature_id_by_accession($dbh, $locus_tag);
		if ($this_feat_id && $this_feat_id ne $feat_id) {
		    print STDERR "ERR: " . $seqo->display_name . " feature at $start..$end (feature_id=$feat_id) has accession ($locus_tag) that is associated with another feature (feature_id=$this_feat_id)\n";
		} elsif (! $this_feat_id) {
		    # need to load accession for feature
		    &insert_feature_accessions($dbh, $feat_id, $locus_tag, $source, "");
		}
	    }
	} elsif (! $feat_id) {
	    my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqo);
	    if (! $feat_id) { die "Failed to load feature" }
	}
	
	&delete_feature_annotations($dbh, $feat_id, {'source' => $source});
	&load_feature_annotations($dbh, $feat_id, $featObj->annotation, $source, 10);
	$locus_tag = "";
	$start = 0;
	$end = 0;
    }
}

print "\n";
exit();

