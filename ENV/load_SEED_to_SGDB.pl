#!/usr/bin/perl -w

## So what this loads is an amalgam of the RAST and SEED downloadable spreadsheets
## Pretty much take the RAST sheet, delete the sequence and aa fields, and
## append the subsystem column.
## Headers: contig_id	feature_id	type	location	start	stop	strand	function	aliases	figfam	evidence_codes	Subsystem

## some of this info, like location, aliases, and evidence codes, is not loaded.
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
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
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
    if ($line =~ /^#/) { warn "Skipping comment line\n$line\n"; next; }
    chomp $line;
    my ($ctg_acc,
	$rast_id,                    
	$feat_type,                  
	$loc_acc,
	$end5,
	$end3,
	$strand,
	$product,   # embedded EC
	$aliases,
	$figfam,
	$ev_codes,
	$subsystem) = split /\t/, $line;

    my $asmbl_id;
    if (defined $ACC->{$ctg_acc}) {
	$asmbl_id = $ACC->{$ctg_acc}->{'asmbl_id'} }
    elsif (grep /^$ctg_acc$/, @current_asmbl_ids) {
	$asmbl_id = $ctg_acc }

    if (! $asmbl_id) { die "Couldn't find asmbl_id for $ctg_acc!!\n"; }

    my @ec;
    while ($product =~ /\(\s*ec[\s\:]+(\d(\.[\d\-]+){0,3})\s*\)/ig) {
	push @ec, $1;
    }

    if ($feat_type eq "peg") {
	$feat_type = "CDS";
    } elsif ($feat_type eq "rna") {
	$feat_type = $product =~ /tRNA/ ? "tRNA" :
	    $product =~ /rRNA/ ? "rRNA" :
	    $product =~ /tmRNA/ ? "tmRNA" :
	    "sRNA";
    } else {
	print "What feat_type for $feat_type/$product? ";
	$feat_type = <STDIN>;
	chomp $feat_type;
    }

    my $feat_name = &make_feat_name($rast_id, $feat_type);
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
    if ($figfam ne "") {
	my %prop = ('feat_name' => $feat_name,
		    'prop_type' => 'FIGfam',
		    'method' => 'RAST',
		    'curated' => 0);
	my @prop_score = ({'score_type' => 'FIGfam',
			   'score' => $figfam});
	&insert_properties($dbh, \%prop, \@prop_score);
    }
}

exit();

sub make_feat_name {
    my $rast_id = shift;
    my $feat_type = shift;

    if ($rast_id =~ /\.\w+\.(\d+)$/) {
	my $index = $1;

	return sprintf "$feat_type%05i", $index;
    } else {
	print "Couldn't find index in $rast_id. Suggest one? \n";
	my $index = <STDIN>;
	chomp $index;
	return sprintf "$feat_type%05d", $index;
    }
}
