#!/usr/bin/perl

use Bio::SeqIO;

my $fileo = Bio::SeqIO->new(-file => $ARGV[0]);

my $seqi;
while (my $seqo = $fileo->next_seq) {
    $seqi++;
    my $locus;
    my @features = $seqo->get_SeqFeatures();
    foreach my $featObj (@features) {
	my @tags = $featObj->get_all_tags;
	if (grep /locus_tag/, @tags) {
	    ($locus) = $featObj->get_tag_values('locus_tag');
	}

	if ($featObj->primary_tag eq "CDS") {
	    my ($end5, $strand) = $featObj->strand > 0 ?
		($featObj->start, "+") : 
		($featObj->end, "-");
	    print "$locus\t$end5\t$strand\t$seqi\n";
	}
    }
}
