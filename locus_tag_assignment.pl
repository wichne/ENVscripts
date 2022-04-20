#!/usr/bin/perl

use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;
use strict;

my %arg;
my $DEBUG = 0;

&getopts('D:u:p:P:i:n:z:s:htS:',\%arg);
if ($arg{h}) { 
    die "USAGE: locus_tag_assignment.pl -D db -p pswd -P locus_tag_prefix [ -i set_id OR -n set_name ] [ -u user -z zerofill_length -s seq_id -t -S source(default=locus_tag)]\n";
}

my $prefix = $arg{'P'};

my $dbh = &connect(\%arg);

my $set_id;
if ($arg{'n'}) {
    $set_id = &set_name_to_id($dbh, $arg{'n'});
} elsif ($arg{'i'}) {
    $set_id = $arg{'i'};
} else { die "Barf;\n"; }
my $zeros = $arg{z} ? $arg{z} : 5;
my $source = $arg{'S'} ? $arg{'S'} : 'locus_tag';
my $seq_aref;
my $i;

if ($arg{'s'}) { 
    $seq_aref = [$arg{'s'}];
    # if we're just assinging to a single sequence, it is likely that in the database
    # there are other locus_tags assigned, so dip into the db and pull out the max 
    # locus_tag and increment from there
    my $max_locus_tag = &get_max_locus_tag($dbh, $set_id);
    if ($max_locus_tag =~ /(\d+)$/) {
	$i = $1;
    } else {
	die "malformed max locus_tag $max_locus_tag for set $set_id\n";
    }
} else {
    $seq_aref = &set_id_to_seq_ids($dbh, $set_id);
}

foreach my $seq_id (@$seq_aref) {
    my $feat_aref = &get_feature_ids_by_seq_id($dbh, $seq_id);
    my $feat_ref = &get_features_by_feature_id($dbh, @$feat_aref);
    #sort feature_ids by position
    foreach my $fid (sort {$feat_ref->{$a}->{'location'}->{$seq_id}->{'feat_min'} <=>
			       $feat_ref->{$b}->{'location'}->{$seq_id}->{'feat_min'} } keys %$feat_ref) {
		if ($feat_ref->{$fid}->{'feat_type'} ne 'CDS' &&
		    $feat_ref->{$fid}->{'feat_type'} !~ /RNA$/) { next }
		$i += 5;
		#assign a locus_tag indexed by 5s.
		my $locus_tag = sprintf $prefix . "_%0${zeros}i", ($i);
		
		# check for existing locus_tag in feature_accessions
		my $existing_locus_ref=&get_accession_by_feature_id($dbh, $fid, $source);
		if (@$existing_locus_ref > 0) {
			foreach my $existing_locus (@$existing_locus_ref) {
				if ($existing_locus eq $locus_tag) { next }
				else {
					warn "$fid: Existing accession '$existing_locus' will be set to '${source}_old' and new accession '$locus_tag' will be inserted as '$source'\n";
					# if it exists rename the source to 'locus_tag_old'
					&set_accession_to_old($dbh, $fid, $source, $existing_locus);
				}
			}
		}
	
		if (!$arg{'t'}) {
			# insert a new row with source 'locus_tag'
			my $locus_tag_i = "INSERT INTO feature_accessions"
			. " (feature_id, source, accession)"
			. " VALUES ($fid, '$source', '$locus_tag')";
			$dbh->do($locus_tag_i) unless ($DEBUG);
		} else {
			print "$seq_id -> $fid -> $locus_tag\n";
		}
    }
}


sub get_max_locus_tag {
    my $dbh = shift;
    my $set_id = shift;
    my $q = "select max(a.accession) from feature_accessions a, seq_feat_mappings m, seq_set_link l"
	. " where source=\"locus_tag\" and m.feature_id=a.feature_id and m.seq_id=l.seq_id"
	. " and l.set_id=$set_id";
    my ($r) = $dbh->selectrow_array($q);
    if ($r) { return $r }
    else { return 0 }
}

sub set_accession_to_old {
    my $dbh = shift;
    my $fid = shift;
    my $source = shift;
    my $accession = shift;

    my $q = "UPDATE feature_accessions set source = '${source}_old'"
	. " WHERE feature_id=$fid AND source='$source'"
	. " AND accession='$accession'";
    my $r = $dbh->do($q);
}
