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
&getopts('D:u:p:i:F:f:s:S:C', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $set_id = $args->{'s'};
my $source = $args->{'S'} or die "Need to provide information source with -s\n";
my $update_coords = $args->{'C'}; 
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

while (my $seqo = $in->next_seq) {

	# if there's no set_id, load a set
	if (! $set_id) {
		warn "Loading a set: $filename\n";
		my $name = $filename;
		$name =~ s/\.[^\.]+$//;
		$set_id = load_sequence_sets($dbh, $name, $seqo->desc);
		if (! $set_id) { die "Can't do a sequence_sets insertion\n"; }
		else { warn "\t$set_id loaded.\n"; }
	}

	# get the seq_id or load the sequence
    my $seq_id = &get_seq_id_by_seq_accession($dbh, $seqo->display_name);
    if (! $seq_id) {
		$seq_id = &load_sequence_SeqFeature($dbh, $seqo);
		if (! $seq_id) { die "Failed loading sequence " . $seqo->display_name . "\n"; }
		else { print STDERR "Sequence " . $seqo->display_id . " loaded: $seq_id\n"; }
		link_seq_to_set($dbh, $set_id, $seq_id);
    }
    print STDERR $seqo->display_name . "|";
	# get the sequence features
    my @features = $seqo->get_SeqFeatures();
    
	# these need to be outside the loop because sometimes gene features have
	# information relevant to CDS features (like locus_tag)
    my ($start, $end, $locus_tag);
    foreach my $featObj (@features) {
		print STDERR ".";
		# Don't load the entire sequence as a feature
		my $thisType = $featObj->primary_tag;
        if ($thisType eq "source") { next }
		if (! grep /^$thisType$/, ("CDS", "rRNA", "tRNA", "sRNA", "misc_RNA") ) { print $thisType; next }
 
		# get the feature coordinates and find any feat_id
		$start = $featObj->start;
		$end = $featObj->end;
		my $feat_id = &get_feature_id_by_coords($dbh, $seq_id, $thisType, 
		                                        $featObj->strand, $start, $end);
		# handle the feature
		# if its a pseudo, load it
		if ($thisType eq "gene") {
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
		
		# Only reach this part if the feature is not a gene

		# This sets the locus_tag to the locus_tag harvested from the gene
		if ($featObj->has_tag('locus_tag')) {
			$locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
		} elsif ($featObj->start eq $start && $featObj->end eq $end && $locus_tag) {
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
				my $this_feat_id = &get_feature_id_by_accession($dbh, $locus_tag, $seq_id);
				if ($this_feat_id && $this_feat_id ne $feat_id) {
					print STDERR "ERR: " . $seqo->display_name 
					. " feature at $start..$end (feature_id=$feat_id) has accession ($locus_tag)"
					. " that is associated with another feature (feature_id=$this_feat_id)\n";
				} elsif (! $this_feat_id) {
					# need to load accession for feature
					&insert_feature_accessions($dbh, $feat_id, $locus_tag, $source, "");
				}
			}
			my ($feat_min, $feat_max) = sort {$a<=>$b} ($start, $end);
			&update_feature_mapping($dbh, $seq_id, $feat_id, {'feat_min'=>$feat_min, 'feat_max'=>$feat_max}) if ($update_coords);
			&delete_feature_annotations($dbh, $feat_id, {'source' => $source});
			&load_feature_annotations($dbh, $feat_id, $featObj, $source, 10);

		} elsif (! $feat_id) {
			$feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqo);
			if (! $feat_id) { die "Failed to load feature" }
			&load_feature_annotations($dbh, $feat_id, $featObj, $source, 10);
		}
		
		$locus_tag = "";
		$start = 0;
		$end = 0;
    }
	print STDERR "\n";
}

exit();

