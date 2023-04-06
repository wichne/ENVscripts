#!/usr/bin/perl -w

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;

my $args = {};
&getopts('f:F:D:p:s:h', $args);

if ($args->{h}) { die "-f file\n-F format\n-D database\n-P password\n-sset_id\n-h [print this]\n"; }

my $filename = $args->{'f'} or die "Need to provide filename with -f\n";
my $database = $args->{'D'} or die "Need to provide database with -D\n";
my $password = $args->{'p'} or die "Need to provide password with -p\n";
my $format = $args->{'F'};
my $set_id = $args->{'s'};

my $in = Bio::SeqIO->new(-file => $filename);

my $dbh = &connect($args);

if (! $set_id) {
    # see if the file prefix is the sequence_sets.name
    if ($filename =~ /(.*)\.[^\.]+$/) {
	my $qname = $1;
	warn("Looking for set_id based on '$qname'...");
        $set_id = set_name_to_id($dbh, $qname);
	if (! $set_id) { warn "No set_id and provided and can't find one based on input filename.\n"; }
	else { warn("\tFound set_id $set_id!");}
    }
}
while (my $seqo = $in->next_seq) {
    my $seq_id;
    if (! $set_id) {
	if ($seqo->display_id && $seqo->description) {
	    my $name = $seqo->display_id;
	    my $desc = $seqo->description;
	    $set_id = load_sequence_sets($dbh, $name, $desc, 0); 
	} else { die "Can't find a proper name/description for this file ($filename) from the header line: "
		     . $seqo->primary_id . ":"
		     . $seqo->display_id . ":"
		     . $seqo->description . "\n"
	}
    }
#    if ($seqo->length >= 2000) {
	$seq_id = load_PrimarySeq($dbh, $seqo);
#    } else { print STDERR $seqo->display_id . " is too short: " . $seqo->length . "\n";
#    }

    if ($set_id) {
	link_seq_to_set($dbh, $set_id, $seq_id);
    }
}

$dbh->disconnect;

exit();
