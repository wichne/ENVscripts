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
&getopts('f:F:D:P:s:h', $args);

if ($args->{h}) { die "-f file\n-F format\n-D database\n-P password\n-sset_id\n-h [print this]\n"; }

my $filename = $args->{'f'} or die "Need to provide filename with -f\n";
my $database = $args->{'D'} or die "Need to provide database with -D\n";
my $password = $args->{'P'} or die "Need to provide password with -P\n";
my $format = $args->{'F'};
my $set_id = $args->{'s'};

my $in = Bio::SeqIO->new(-file => $filename, -format => 'swiss');

my $dbh = &connect($args);

# if (! $set_id) {
#     # see if the file prefix is the sequence_sets.name
#     if ($filename =~ /(.*)\.[^\.]+$/) {
# 	my $qname = $1;
# 	warn("Looking for set_id based on '$qname'...");
#         $set_id = set_name_to_id($dbh, $qname);
# 	if (! $set_id) { die "No set_id and provided and can't find one based on input filename.\n"; } else { warn("\tFound set_id $set_id!");}
#     }
# }
while (my $seqo = $in->next_seq) {
    my $seq_id;
    my $annotobj = $seqo->annotation;
    # this gives us
    # comment
    # secondary_accession
    # date_changed
    # keyword
    # evidence
    # dblink
    # reference
    # gene_name
    # seq_update
    my %ANN;
    my ($dblink) = $annotobj->get_Annotations('dblink');
    my (@link) = split(/\n/, $dblink);
    foreach my $l (@link) {
	if ($l =~ /link to (\S+) .* database (eggNOG)/) { push @{$ANN{$2}}, $1 }
	if ($1 =~ /link to (\S+) .* database (GO)/) { push @{$ANN{$2}}, $1 }
	if ($1 =~ /link to (\S+) .* database (InterPro)/) { push @{$ANN{$2}}, $1 }
	if ($1 =~ /link to (\S+) .* database (Pfam)/) { push @{$ANN{$2}}, $1 }
	if ($1 =~ /link to (\S+) .* database (SUPFAM)/) { push @{$ANN{$2}}, $1 }
	if ($1 =~ /link to (\S+) .* database (TIGRfams)/) { push @{$ANN{$2}}, $1 }
    }

    my @gene_names = $annotobj->get_Annotations('gene_name');
    foreach my $gn (@gene_names) {
		my @vals = $val->get_all_values;
		print join("\n+\t", @vals), "\n";
	    } else {
		print $val->as_text() . "\n";
	    }
	}
    }
#    if ($seqo->length >= 2000) {
#	$seq_id = &load_PrimarySeq($dbh, $seqo);
#    } else { print STDERR $seqo->display_id . " is too short: " . $seqo->length . "\n";
#    }

#    if ($set_id) {
#	&link_seq_to_set($dbh, $set_id, $seq_id);
#    }
}

$dbh->disconnect;

exit();
