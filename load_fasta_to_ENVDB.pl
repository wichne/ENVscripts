#!/usr/bin/perl -w

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

my $args = {};
&getopts('f:F:D:P:s:h', $args);

if ($args->{h}) { die "-f file\n-F format\n-D database\n-P password\n-sset_id\n-h [print this]\n"; }

my $filename = $args->{'f'} or die "Need to provide filename with -f\n";
my $database = $args->{'D'} or die "Need to provide database with -D\n";
my $password = $args->{'P'} or die "Need to provide password with -P\n";
my $format = $args->{'F'};
my $set_id = $args->{'s'};

my $in = Bio::SeqIO->new(-file => $filename);

my $dbh = &connect($args);

while (my $seqo = $in->next_seq) {
    my $seq_id;
#    if ($seqo->length >= 2000) {
	$seq_id = &load_PrimarySeq($dbh, $seqo);
#    } else { print STDERR $seqo->display_id . " is too short: " . $seqo->length . "\n";
#    }

    if ($set_id) {
	&link_seq_to_set($dbh, $set_id, $seq_id);
    }
}

$dbh->disconnect;

exit();
