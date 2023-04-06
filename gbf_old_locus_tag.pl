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

my $args = {};
&getopts('D:u:p:i:F:s:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $source = $args->{'s'} ? $args->{'s'} : "PNNL";

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

    foreach my $featObj (@features) {
	if ($featObj->primary_tag eq "gene") {
	    if ($featObj->has_tag('locus_tag')) {
		my $locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
		if ($featObj->has_tag('old_locus_tag')) {
		    my $other_acc = [$featObj->get_tag_values('old_locus_tag')]->[0];
		    my $feat_id = &get_feature_id_by_accession($dbh, $locus_tag);
		    if (! $feat_id) { print STDERR "Why no feature with locus_tag '$locus_tag' ?\n"; }
		    else {
			my $iok = &insert_feature_accessions($dbh, $feat_id, $other_acc, $source);
			if (! $iok) { warn "Insert of $other_acc for $feat_id failed.\n"; }
		    }
		}
	    }
	}
    }
}

    print "\n";
    exit();

