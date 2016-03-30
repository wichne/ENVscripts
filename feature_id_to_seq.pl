#!/usr/bin/perl

#Edited by Jen 05/08/2015 to change sequence identifier to database,feature id and organism 

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

&getopts('D:i:o:');
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'access');

my $fid = $opt_i;
my $outfh;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
}

our $SEQ;

while (my $fid = <STDIN>) {
    chomp $fid;
    my $protref = &get_features_by_feature_id($dbh, $fid);
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
    my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
    my $db = $opt_D;
    my $acc = "$db|$fid|gnl|" . $locus->[0];
    
    my @seq_ids = keys (%{$protref->{$fid}->{'location'}});
    if (@seq_ids > 1) {
	warn "This feature $fid maps to multiple seq_ids @seq_ids";
    }
	
    my ($seq_id, $min, $max, $strand) = ($protref->{$fid}->{'location'}->{$seq_ids[0]}->{'seq_id'},
					 $protref->{$fid}->{'location'}->{$seq_ids[0]}->{'feat_min'},
					 $protref->{$fid}->{'location'}->{$seq_ids[0]}->{'feat_max'},
					 $protref->{$fid}->{'location'}->{$seq_ids[0]}->{'strand'});
    if (!($min && $max) || ($min >= $max)) { warn "Why $min/$max for $fid, $acc?\n";}
    my $seq = &get_subseq($seq_id, $min, $max, $strand);
    
#     if (length($protref->{$fid}->{'product'}) == 0 ||
# 	!defined($protref->{$fid}->{'product'})) {
# 	warn "No product string for feature $fid ($acc). Skipping...";
# 	next;
#     }
    my $setdesc = &get_set_desc($dbh, $protref->{$fid}->{'set'}->{'id'});
    my $desc = "[$setdesc]";
    while (my ($rank,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	my ($source) = sort {$a<=>$b} keys %$vref;
	foreach my $qual (keys %{$vref->{$source}}) {
	    $desc .= " $qual=$vref->{$source}->{$qual}->[0];";
	}	
    }

    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $seq);
    $outfh->write_seq($seqo);
}

sub get_set_desc {
    my $dbh = shift;
    my $setid = shift;
    my $ssq = "SELECT description FROM sequence_sets WHERE set_id=$setid";
    return $dbh->selectcol_arrayref($ssq)->[0];
}

sub get_mainroles {
    my $dbh = shift;
    my $mrq = "select * from egad.mainrole";
    return $dbh->selectall_hashref($mrq, 'id');
}

sub get_subroles {
    my $dbh = shift;
    my $mrq = "select * from egad.subrole";
    return $dbh->selectall_hashref($mrq, 'id');
}

sub get_subseq {
    my ($seq_id, $min, $max, $strand) = @_;
    $SEQ = &get_sequence_by_seq_id($dbh, $seq_id);
    my $subseq = $SEQ->{$seq_id}->trunc($min, $max);
    if ($strand < 0 || $strand eq "-") {
	return $subseq->revcom->seq;
    } else {
	return $subseq->seq;
    }
}
