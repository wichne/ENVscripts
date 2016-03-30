#!/usr/bin/perl
use Bio::SeqIO;
use Getopt::Std;
use lib $ENV{SCRIPTS};
use ENV;

my %arg;
&getopts('D:p:u:f:', \%arg);

my $dbh = &connect(\%arg);

my $file = $arg{'f'};
my $fileo = Bio::SeqIO->new(-file => $file);

my $product_q = "update sequence_features set product=? where feature_id=?";
my $sth1 = $dbh->prepare($product_q);

my $pseudo_q = "update seq_feat_mappings m, seq_set_link l, sequence_sets s"
    . " SET m.pseudo=5 where m.feature_id=?"
    . " AND m.seq_id=l.seq_id and s.set_id=l.set_id and s.is_current=1";
my $sth2 = $dbh->prepare($pseudo_q);

while (my $rec = $fileo->next_seq) {
    my $acc = $rec->display_name;
    my $fid = &get_feature_id_by_accession($dbh, $rec->display_id);
    if (! $fid) { warn "No feature_id for $acc\n"; next; }
    my $pep = $rec->seq;

    $sth1->execute($pep, $fid);
#    $sth2->execute($fid);
}


