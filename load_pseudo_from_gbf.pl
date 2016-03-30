#!/usr/bin/perl
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Std;
use lib $ENV{SCRIPTS};
use ENV;
use strict;

my %arg;
&getopts('D:p:u:g:f:h', \%arg);
if ($arg{'h'} || ! $arg{'p'}) {
    die "USAGE: load_pseudo_from_gbf -D db -p pswd [-u user] -g gbf_file -f aa_fasta_file";
}
my $dbh = &connect(\%arg);

# crack open the Genbank file:
my $file = $arg{'g'};
my $fileo = Bio::SeqIO->new(-file => $file, -format => 'genbank');

# crack open the peptide fasta file
my $pep = $arg{'f'};
my $dbo = Bio::DB::Fasta->new($pep);

while (my $seqo = $fileo->next_seq) {
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence into assembly, stan and asmbl_data
    my $seqObj = shift @features;
    my $seq_id = get_seq_id_by_seq_accession($dbh, $seqo->display_name);
    if (! $seq_id ) { die "Couldn't get seq_id from " . $seqo->display_name . "\n"; }

    my ($start, $end, $locus_tag);
    foreach my $featObj (@features) {
	if ($featObj->primary_tag eq "gene") {
	    if ($featObj->has_tag('pseudo')) {
		my $fid = &get_feature_id_by_accession($dbh, $featObj->get_tag_values('locus_tag'));
		my $current = is_current($dbh, $fid);
		my ($locus_tag) = $featObj->get_tag_values('locus_tag');
		my $aa = $dbo->get_Seq_by_id($locus_tag)->seq;
		if ($fid && $current) {
		    print "$locus_tag already loaded\n";
		    # update sequence_features.product field
		    # set the feature to pseudo
		    my $q = "update seq_feat_mappings set pseudo=5 where seq_id=$seq_id and feature_id=$fid";
		    $dbh->do($q);
		    next;
		} elsif ($fid) {
		    print "Feature not on current molecule. Loading seq_feat_mappings and updating product...\n";
		    load_seq_feat_mappings($dbh, $fid, $seq_id, $featObj);
		    &update_product($dbh, $fid, $aa);
		    my $q = "update seq_feat_mappings set pseudo=5 where seq_id=$seq_id and feature_id=$fid";
		    $dbh->do($q);
		} else {
		    $seqObj->add_tag_value('translation', $aa);
		    $fid = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj, undef, "JGI");
		    print "inserted ", $featObj->display_name, " as $fid $locus_tag at " . $featObj->start . " / " . $featObj->end . "\n";
		}
	    } else {
#		my $q = "update seq_feat_mappings set pseudo=0 where seq_id=$seq_id and feature_id=$fid";
#		$dbh->do($q);
	    }
	} 
    }
}
