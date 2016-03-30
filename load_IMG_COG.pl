#!/usr/bin/perl

# Edited by JMo on 12/09/2014 to accomodate the COG data grepped from JGI_info file (*.info.xls)

use Getopt::Std;
use DBI;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

my %arg;
&getopts('D:u:p:i:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
if (! $infile) { die "No COG file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}

open my $in, $infile or die "Can't open $infile: $!\n";

while (my $line = <$in>) {
    chomp $line;
#Commented out from a previous format which is unknown
   #my ($acc, $len, $lo, $hi, $score, $cog) = split/\,/, $line;
   #next if (! $score);
   # my @acc = split/\|/, $acc;
   # my $fid = &get_feature_id_by_accession($dbh, $acc[0]);
   # my $len = &get_feature_product_length($dbh, $fid);

    my ($jgi, $acc, $COG_acc, $desc, $evalue) = split/\t/, $line;
    my $fid = &get_feature_id_by_accession($dbh, $acc);
    if (!$fid) {die "Why no fid for $acc?";}
    my $len = &get_feature_product_length($dbh, $fid);
    my $ev_i = "INSERT INTO feature_evidence (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
	. " VALUES (?, '1', ?, 'IMG_pipeline', 'COG', ?, '0', '0', '0', ?)";
    my $evstg = $dbh->prepare($ev_i);
    $evstg->execute($fid, $len, $COG_acc, $evalue);
}

sub get_feature_product_length {
    my $dbh = shift;
    my $fid = shift;
    my $fq = "SELECT length(product) FROM sequence_features WHERE feature_id=$fid";
    my @l = $dbh->selectrow_array($fq);
    if ($l[0]==0) {
	my $mq = "SELECT (feat_max - feat_min + 1)/3 from seq_feat_mappings m where feature_id=$fid";
	@l = $dbh->selectrow_array($mq);
    }
    return $l[0];
}
