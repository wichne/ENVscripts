#!/usr/bin/perl -w

use DBI;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{SCRIPTS};
use SGDB;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:a:', $args);

my $filename = $args->{'i'};
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $asmbl_id = $args->{'a'};
my $host = "bladerunner:3306";
$dbh = DBI->connect("dbi:mysql:$db:$host",
		       $user, $pswd, {'AutoCommit'=> 0});
my $in;
if ($filename) {
    open $in, $filename or die "Can't open '$filename': $!\n";
} else {
    print STDERR "Expecting input from STDIN...\n";
    $in = *STDIN;
}

my ($ACC, @current_asmbl_ids);
if (! $asmbl_id) {
    my $acc_q = "SELECT s.asmbl_id, a.accession"
	. " FROM stan s, asmbl_data d, asmbl_accessions a"
	. " WHERE s.iscurrent =1"
	. " AND d.id = s.asmbl_data_id"
	. " AND a.asmbl_data_id = d.id";
    $ACC = $dbh->selectall_hashref($acc_q, 'accession'); # hash ref

    foreach my $ref(values %$ACC) {
	push @current_asmbl_ids, $ref->{'asmbl_id'};
    }
}

while (my $line = <$in>) {
    if ($line =~ /^Feature ID/) { next } # header line
    if ($line !~ /^fig/) { warn "Line did not start with RAST prefix 'fig'\n$line\n"; next; }
    chomp $line;
    my (
	$rast_id,
	$feat_type,
	$ctg_acc,
	$end5,
	$end3,
	$frame,
	$strand,
	$length,
	$product,   # embedded EC
	$subsystem,
	$ncbi_gi_locus) = split /\t/, $line;

    my $asmbl_id;
    if (defined $ACC->{$ctg_acc}) {
	$asmbl_id = $ACC->{$ctg_acc}->{'asmbl_id'} }
    elsif (grep /^$ctg_acc$/, @current_asmbl_ids) {
	$asmbl_id = $ctg_acc }

    my @ec;
    while ($product =~ /\(\s*ec[\s\:]+(\d(\.[\d\-]+){0,3})\s*\)/ig) {
	push @ec, $1;
    }

    my $feat_name = &get_next_feat_name($dbh, {'feat_type' => $feat_type} );
    my %asm_feature = ('asmbl_id' => $asmbl_id,
		       'feat_type' => $feat_type,
		       'feat_name' => $feat_name,
		       'end5' => $end5,
		       'end3' => $end3,
		       );
    &insert_feature($dbh, \%asm_feature);

    my %ident = ('feat_name' => $feat_name,
		 'product' => $product,
		 'ec' => join(" ", @ec),
#		 'comment' => $comment
	);
    &insert_ident($dbh, \%ident);

    # insert prot acc
    my ($prefix, $acc) = split /\|/, $rast_id;
    if (! $acc) { $acc = $prefix; $prefix = "" }
    my %feat_acc = ('feat_name' => $feat_name,
		    'source' => 'RAST',
		    'prefix' => $prefix,
		    'accession' => $acc);
    &insert_feature_accessions($dbh, \%feat_acc);

    # insert subsystem into properties
    if ($subsystem ne "- none -") {
	my %prop = ('feat_name' => $feat_name,
		    'prop_type' => 'SEED subsystem',
		    'method' => 'RAST',
		    'curated' => 0);
	my @prop_score = ({'score_type' => 'subsystem',
			   'score' => $subsystem});
	&insert_properties($dbh, \%prop, \@prop_score);
    }
}

exit();
