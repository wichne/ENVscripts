#!/usr/bin/perl -w

use DBI;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib "/home/nels329/devel";
use SGDB;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:a:', $args);

my $filename = $args->{'i'};
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $asmbl_id = $args->{'a'};

$dbh = DBI->connect("dbi:mysql:server=localhost;database=$db",
		       $user, $pswd, {'AutoCommit'=> 0});
my $in;
if ($filename) {
    open $in, $filename or die "Can't open '$filename': $!\n";
} else {
    print STDERR "Expecting input from STDIN...\n";
    $in = *STDIN;
}

if (! $asmbl_id) {
    my $current_asmbl_ids = &current_asmbl_ids($dbh, $db);
    if (@$current_asmbl_ids == 1) {
	$asmbl_id = $current_asmbl_ids->[0]->[0];
    } else {
	print STDERR "Please specify asmbl_id with -a. Current asmbl_ids:\n";
	print STDERR "$_->[0]\n" foreach (@$current_asmbl_ids);
	die;
    }
}

while (my $line = <$in>) {
    chomp $line;
    my ($end5,
	$end3,
	$feat_type,
	$locus,
	$com_name,
	$gene_sym,
	$ec,
	$comment) = split /\t/, $line;
    my $feat_name = &get_next_feat_name($dbh, {'feat_type' => $feat_type} );
    my %asm_feature = ('asmbl_id' => $asmbl_id,
		       'feat_type' => $feat_type,
		       'feat_name' => $feat_name,
		       'end5' => $end5,
		       'end3' => $end3,
		       );
    &insert_feature($dbh, \%asm_feature);

    my %ident = ('feat_name' => $feat_name,
		 'locus' => $locus,
		 'product' => $com_name,
		 'gene_sym' => $gene_sym,
		 'ec' => $ec,
		 'comment' => $comment);
    &insert_ident($dbh, \%ident);

}

exit();
