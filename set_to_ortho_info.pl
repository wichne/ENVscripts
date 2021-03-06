#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

&getopts('D:i:n:o:');
my $db = $opt_D;
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'access');

my $setid = $opt_i;
my $setname = $opt_n;
my $outfh;
if ($opt_o) {
    open($outfh, ">$opt_o");
} else {
    $outfh = \*STDOUT;
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $featref = &get_seq_features_by_set_id($dbh, $setid, "CDS");

foreach my $fid (keys %$featref) {
#    my $acc = join "|", values %{$featref->{$fid}->{'accessions'}};
#    my $locus_tag = get_locus_tag_by_feature_id($dbh, $fid);
    my $locus_tag = $featref->{$fid}->{'accessions'}->{'locus_tag'};
#    my $db = $opt{D};
    my $acc = "$db|$fid|$locus_tag";
#    my $acc = join "|", ($db, $fid, $locus_tag->[0]);
    my $dir = $featref->{$fid}->{'location'}->{'strand'} > 0 ? "+" : "-";
    print $outfh join("\t", ($acc,
			     $featref->{$fid}->{'location'}->{'feat_min'},
			     $dir,
			     $featref->{$fid}->{'location'}->{'seq_id'}))
	. "\n";
}
close $outfh;
exit();
