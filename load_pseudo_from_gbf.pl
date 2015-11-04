#!/usr/bin/perl
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Std;
use lib $ENV{SCRIPTS};
use ENV;

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
		my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj, undef, "JGI");
		my $aa = $dbo->get_Seq_by_id($featObj->display_id);
		print "inserted $feat_id " . $featObj->get_tag_values('locus_tag') . " at " . $featObj->start . " / " . $featObj->end . "\n";
	    }
	}
    }
}
