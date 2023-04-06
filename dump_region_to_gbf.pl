#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocationI;

my $opt = {};
&getopts('D:i:o:c:C:', $opt);
my $db = $opt->{D};
my $seqid = $opt->{i};
my $outfile = $opt->{o};
my $locoord = int($opt->{c});
my $hicoord = int($opt->{C});

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'mySQL@cce55');

my $OUT;
if ($outfile) {
    my $o = $outfile;
    $OUT = Bio::SeqIO->new(-file => ">$o",
			   -format => 'genbank');
}
if (! $outfile) {
    $OUT = Bio::SeqIO->new(-fh => \*STDOUT,
			   -format => 'genbank');
}

my $SeqObj = &seq_id_to_SeqObj($dbh, $seqid);
add_features_to_SeqObj($dbh, $seqid, $SeqObj);
# get the subsequence and make a new seq obj
my $subseq = $SeqObj->subseq($locoord, $hicoord);
my $ssobj = Bio::Seq->new(-display_id => $SeqObj->display_id . "_${locoord}_${hicoord}", -seq => $subseq);
my @features = $SeqObj->get_SeqFeatures();
foreach my $featObj (@features) {
    my $start = $featObj->start;
	my $end = $featObj->end;
    if (($locoord < $end && $hicoord > $end) 
        || ($locoord < $start && $hicoord > $start) 
        || ($start < $locoord && $end > $hicoord)) {
        my $thisFeat = $featObj;
        $thisFeat->start($featObj->start - ($locoord-1));
        if ($thisFeat->start < 0) { 
            $thisFeat->start('<1');
        }
        if ($featObj->end > $hicoord) {
            $thisFeat->end('>' . ($hicoord - ($locoord - 1)));
        } else {
            $thisFeat->end($featObj->end - ($locoord-1));
        }
        if ($thisFeat->primary_tag eq "source") {
            $thisFeat->start(1);
            $thisFeat->end($hicoord - ($locoord -1));
        }
        $ssobj->add_SeqFeature($thisFeat);
    }
}

$OUT->write_seq($ssobj);

