#!/usr/bin/perl
#Edited 05/17/2015 by Jen in order to incorporate only the databasename and feature id and organism description
#in sequence headers

use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);
use strict;

my %opt;
&getopts('D:i:n:o:A:ah', \%opt);

if ($opt{h}) {
    die "USAGE\nset_to_pep.pl -D db [ -i set_id -n set_name ] [-o output_file -A source_for_accession -a (include annotation) ]\n";
}

my $db = $opt{D};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'mySQL@cce55');

my $setid = $opt{i};
my $setname = $opt{n};
my $outfh;
my $pseudo_file;
if ($opt{o}) {
    $outfh = Bio::SeqIO->new(-file => ">$opt{o}",
			     -format => 'fasta');
    $pseudo_file = $opt{o} . ".pseudo";
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
    $pseudo_file = $$ . ".pseudo";
}
open (my $PSEUDO, ">$pseudo_file");

#my $setref = &get_seq_sets($dbh);
## Make sure we have a set id and a set name
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

## get the features from the setid
my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");

## get the seqids from the setid
my $seqids = &set_id_to_seq_ids($dbh, $setid);

## get the set description
my $setdesc = &get_set_desc($dbh,$setid);

## get the sequences by the seqids
our $SEQ = &get_sequence_by_seq_id($dbh, @$seqids);

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
    my $acc;
    if ($opt{A}) {
	if ($opt{A} eq "all") {
	    $acc = "gnl|${db}_$fid";
	    foreach my $src (keys %{$protref->{$fid}->{'accessions'}}) {
		$acc .= "|" . $protref->{$fid}->{'accessions'}->{$src};
	    }
	} elsif (defined $protref->{$fid}->{'accessions'}->{$opt{A}}) {
	    $acc = $protref->{$fid}->{'accessions'}->{$opt{A}};
	} else {
	    $acc = "gnl|${db}_$fid";
	}
    } else {
	#my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
	$acc = "${db}_$fid";
    }

    my $gene_length = $protref->{$fid}->{'location'}->{'feat_max'} - $protref->{$fid}->{'location'}->{'feat_min'} + 1;
    my $calc_prot_length = $gene_length/3 - 1; # gene length should include stop codon but protein should not
    my $prot_length = length($protref->{$fid}->{'product'});
    if (($protref->{$fid}->{'pseudo'} == 0 &&
	 abs($prot_length - $calc_prot_length) > 0) ||
	($protref->{$fid}->{'pseudo'} > 0 &&
	abs($prot_length - $calc_prot_length) > 2)) {
	warn "MISMATCHED TRANSLATION LENGTH: $fid $acc gene_length: $gene_length, calc_prot_length: $calc_prot_length, prot_length: $prot_length\n";
    }
    if (length($protref->{$fid}->{'product'}) == 0 ||
	!defined($protref->{$fid}->{'product'})) {# ||
#	$protref->{$fid}->{'product'} =~ /\*./)) {
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
				 -desc => 'temp',
				 -seq => $seq);
	my $complete = $min_partial || $max_partial ? 0 : 1; 
	my $protobj = $seqo->translate(-complete => $complete,
				       -frame    => $phase,
				       -codontable_id => 11);
	$protref->{$fid}->{'product'} = $protobj->seq;
    }
    my $desc = "[$setdesc]";
    if ($opt{a}) {
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    $desc .= " $qual=$vref->{$k}->[0]->{value};";
	}
	if ($protref->{$fid}->{'location'}->{'pseudo'} > 0) {
	    $desc .= " pseudogene=" . $ENV::PSEUDOKEY{$protref->{$fid}->{location}->{pseudo}} . ";";
	}
    }
    if ($protref->{$fid}->{'product'} =~ /\*./) {
	warn "INTERNAL STOP: $fid $acc $protref->{$fid}->{pseudo}\n";
	printf $PSEUDO "$fid\t%s\t%s\t%s\n", ($protref->{$fid}->{'location'}->{'seq_id'},
					      $protref->{$fid}->{'location'}->{'feat_min'},
					      $protref->{$fid}->{'location'}->{'feat_max'});
#	next;
    }
    if ($protref->{$fid}->{'product'} =~ /\*$/) {
	warn "TRIM TERMINAL STOP: $fid $acc\n";
    }
    if ($protref->{$fid}->{pseudo}) {
	$desc .= " PSEUDO";
    }
    if ($protref->{$fid}->{'location'}->{'min_partial'} ||
	$protref->{$fid}->{'location'}->{'max_partial'}) {
	$desc .= " PARTIAL";
    }

    $protref->{$fid}->{'product'} =~ s/\s//g;
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
