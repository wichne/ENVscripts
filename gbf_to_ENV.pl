#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib "/scripts/";
use ENV;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:f:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};

# crack open the Genbank file:
my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

# crack open the peptide fasta file
my $pep = $arg{'f'};
my $dbo = Bio::DB::Fasta->new($pep);

$dbh = DBI->connect("dbi:mysql:host=$host;database=$db",
		       $user, $pswd, {'AutoCommit'=> 0});
while (my $seqo = $in->next_seq) {
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence into assembly, stan and asmbl_data
    my $seqObj = shift @features;
    my $seq_id = &load_sequence_SeqFeature($dbh, $seqObj);
    print STDERR "Sequence " . $seqObj->display_id . " loaded: $seq_id\n";
    my ($start, $end, $locus_tag);
    foreach my $featObj (@features) {
	if ($featObj->primary_tag eq "gene") {
	    if ($featObj->has_tag('locus_tag')) {
		$locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
		$start = $featObj->start;
		$end = $featObj->end;
	    }
	    if ($featObj->has_tag('pseudo')) {
#		my $feat_id = &get_feature_id_by_accession($dbh, $featObj->get_tag_values('locus_tag'));
		my $aa = $dbo->get_Seq_by_id($featObj->display_name);
#		if ($feat_id) {
#		    print "$locus_tag already loaded\n";
#		    # update sequence_features.product field
#		    &update_product($dbh, $feat_id, $aa);
#		} else {
		    $seqObj->add_tag_value('translation', $aa);
		    $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj, undef, "JGI");
#		    print "inserted $feat_id " . $featObj->get_tag_values('locus_tag') . " at " . $featObj->start . " / " . $featObj->end . "\n";
#		}
	    }
	    next;
	}
	
	if ($featObj->start eq $start && $featObj->end eq $end) {
	    $featObj->set_attributes(-tag => {'locus_tag' => $locus_tag });
	}
	my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj);
	&load_feature_annotations($dbh, $feat_id, $featObj->annotation);
    }
}

print "\n";
exit();

