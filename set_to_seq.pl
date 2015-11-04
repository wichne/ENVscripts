#!/usr/bin/perl

#Edited by Jen 05/08/2015 to change sequence identifier to database,feature id and organism 

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

&getopts('D:i:n:o:A:ah');

if ($opt_h) {
    die "USAGE\nset_to_seq.pl -D db [ -i set_id -n set_name ] [-o output_file -A source_for_accession -a (include annotation) ]\n";
}

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'access');

my $setid = $opt_i;
my $setname = $opt_n;
my $outfh;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

if (! $setid) { die "Can't get setid\n"; }

my $protref = &get_seq_features_by_set_id($dbh, $setid);
my $seqids = &set_id_to_seq_ids($dbh, $setid);
our $SEQ = &get_sequence_by_seq_id($dbh, @$seqids);

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {

    my $acc;
    if ($opt_A) {
	if (defined $protref->{$fid}->{'accessions'}->{$opt_A}) {
	    $acc = $protref->{$fid}->{'accessions'}->{$opt_A};
	    $acc =~ s/.*\|//;
	} else {
	    $acc = $fid;
	}
    } else {
	my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
	my $db = $opt_D;
	$acc = "$db|$fid|$locus";
    }
    

    my ($seq_id, $min, $max, $strand) = ($protref->{$fid}->{'location'}->{'seq_id'},
					 $protref->{$fid}->{'location'}->{'feat_min'},
					 $protref->{$fid}->{'location'}->{'feat_max'},
					 $protref->{$fid}->{'location'}->{'strand'});
    if (!($min && $max) || ($min >= $max)) { warn "Why $min/$max for $fid, $acc?\n";}
    my $seq = &get_subseq($seq_id, $min, $max, $strand);
    
#     if (length($protref->{$fid}->{'product'}) == 0 ||
# 	!defined($protref->{$fid}->{'product'})) {
# 	warn "No product string for feature $fid ($acc). Skipping...";
# 	next;
#     }
    my $setdesc = &get_set_desc($dbh, $setid);
    my $desc = "[$setdesc]";
    if ($opt_a) {
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    $desc .= " $qual=$vref->{$k}->[0]->{value};";
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
    my $subseq = $SEQ->{$seq_id}->trunc($min, $max);
    if ($strand < 0 || $strand eq "-") {
	return $subseq->revcom->seq;
    } else {
	return $subseq->seq;
    }
}
