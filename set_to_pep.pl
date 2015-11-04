#!/usr/bin/perl
#Edited 05/17/2015 by Jen in order to incorporate only the databasename and feature id and organism description
#in sequence headers

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

&getopts('D:i:n:o:A:ah');

if ($opt_h) {
    die "USAGE\nset_to_pep.pl -D db [ -i set_id -n set_name ] [-o output_file -A source_for_accession -a (include annotation) ]\n";
}

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'access');

my $setid = $opt_i;
my $setname = $opt_n;
my $outfh;
my $pseudo_file;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
    $pseudo_file = $opt_o . ".pseudo";
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
    $pseudo_file = $$ . ".pseudo";
}
open (my $PSEUDO, ">$pseudo_file");

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");
my $seqids = &set_id_to_seq_ids($dbh, $setid);
our $SEQ = &get_sequence_by_seq_id($dbh, @$seqids);

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
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

    if (length($protref->{$fid}->{'product'}) == 0 ||
	!defined($protref->{$fid}->{'product'} ||
	$protref->{$fid}->{'product'} =~ /\*./)) {
#	warn "No product string for feature $fid ($acc).";
	my ($seq_id, $min, $max, $strand,
	    $min_partial, $max_partial, $phase) = ($protref->{$fid}->{'location'}->{'seq_id'},
						   $protref->{$fid}->{'location'}->{'feat_min'},
						   $protref->{$fid}->{'location'}->{'feat_max'},
						   $protref->{$fid}->{'location'}->{'strand'},
						   $protref->{$fid}->{'location'}->{'min_partial'},
						   $protref->{$fid}->{'location'}->{'max_partial'},
						   $protref->{$fid}->{'location'}->{'phase'});
	if (!($min && $max) || ($min >= $max)) { warn "Why $min/$max for $fid, $acc?\n";}
	my $seq = &get_subseq($seq_id, $min, $max, $strand);
	my $seqo = Bio::Seq->new(-display_id => $acc,
				 -desc => $desc,
				 -seq => $seq);
	my $complete = $min_partial || $max_partial ? 0 : 1; 
	my $protobj = $seqo->translate(-complete => $complete,
				       -frame    => $phase,
				       -codontable_id => 11);
	if ($protobj->seq !~ /\*./) {
	    $protref->{$fid}->{'product'} = $protobj->seq;
	} else {
	    print $PSEUDO "$fid\t$seq_id\t$min\t$max\t" . $SEQ->{$seq_id}->length . "\n";
	    next;
	}
    }
    my $setdesc = &get_set_desc($dbh,$setid);
    my $desc = "[$setdesc]";
    if ($opt_a) {
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    $desc .= " $qual=$vref->{$k}->[0]->{value};";
	}
    }

    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $protref->{$fid}->{'product'});
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
