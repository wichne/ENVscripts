#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib "/share/scripts/";
use ENV;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};

my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
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

