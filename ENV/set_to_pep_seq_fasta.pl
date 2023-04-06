#!/usr/bin/perl
#Edited 05/17/2015 by Jen in order to incorporate only the databasename and feature id and organism description
#in sequence headers

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

my %opt;
&getopts('D:i:n:o:sA:ah', \%opt);

if ($opt{h}) {
    die "USAGE\nset_to_pep.pl -D db [ -i set_id -n set_name ] [-o output_file -A source_for_accession -a (include annotation) ]\n";
}

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt{D}", 'access', 'access');

my $setid = $opt{i};
my $setname = $opt{n};
my $db = $opt{D};

if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }
if ($setid && ! $setname) { $setname = &set_id_to_set_name($dbh, $setid) }
my $pep_file = sprintf "%s.%s.pep", ($db, $setname);
my $seq_file = sprintf "%s.%s.seq", ($db, $setname);
my $pseudo_file = sprintf "%s.%s.pseudo", ($db, $setname);
my $pepfh = Bio::SeqIO->new(-file => ">$pep_file",
			 -format => 'fasta');
my $seqfh = Bio::SeqIO->new(-file => ">$seq_file",
			 -format => 'fasta');
open (my $PSEUDO, ">$pseudo_file");


my $protref = &get_seq_features_by_set_id($dbh, $setid);
my $seqids = &set_id_to_seq_ids($dbh, $setid);
my $setdesc = &get_set_desc($dbh,$setid);
our $SEQ = &get_sequence_by_seq_id($dbh, @$seqids);

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
    my $acc;
    if ($opt{A}) {
	if (defined $protref->{$fid}->{'accessions'}->{$opt{A}}) {
	    $acc = $protref->{$fid}->{'accessions'}->{$opt{A}};
	    $acc =~ s/.*\|//;
	} else {
	    $acc = $fid;
	}
    } else {
	my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
	my $db = $opt{D};
	$acc = "$db|$fid|$locus";
    }

    ###  Get the gene sequence  ###
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
    my $gene_length = $protref->{$fid}->{'location'}->{'feat_max'} - $protref->{$fid}->{'location'}->{'feat_min'} + 1;
    my $gene_length2 = length($seq);
    if ($gene_length != $gene_length2) { die "Do math right!! Calculated length is not equal to actual length for $fid\n" }

    ### Set up the description  ###
    my $desc = "[$setdesc]";
    if ($opt{a}) {
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    $desc .= " $qual=$vref->{$k}->[0]->{value};";
	}
	if ($protref->{$fid}->{'location'}->{'pseudo'} > 0) {
	    $desc .= " pseudogene=" . $PSEUDOKEY{$protref->{$fid}->{location}->{pseudo}} . ";";
	}
    }
    if ($protref->{$fid}->{pseudo}) {
	$desc .= " PSEUDO";
    }
    if ($protref->{$fid}->{'location'}->{'min_partial'} ||
	$protref->{$fid}->{'location'}->{'max_partial'}) {
	$desc .= " PARTIAL";
    }
    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $seq);

    ###  Get the protein sequence  ###
    if ($protref->{$fid}->{'feat_type'} eq "CDS") {
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
	    my $complete = $min_partial || $max_partial ? 0 : 1; 
	    my $protobj = $seqo->translate(-complete => $complete,
					   -frame    => $phase,
					   -codontable_id => 11);
	    $protref->{$fid}->{'product'} = $protobj->seq;
	}
	if ($protref->{$fid}->{'product'} =~ /\*./) {
	    warn "INTERNAL STOP: $fid $acc $protref->{$fid}->{pseudo}\n";
	    print $PSEUDO "$fid\t$seq_id\t$min\t$max\n";
#	next;
	}
	if ($protref->{$fid}->{'product'} =~ /\*$/) {
	    warn "TRIM TERMINAL STOP: $fid $acc\n";
	}
	
	###  Print out sequences  ###
	$protref->{$fid}->{'product'} =~ s/\s//g;
	my $pepo = Bio::Seq->new(-display_id => $acc,
				 -desc => $desc,
				 -seq => $protref->{$fid}->{'product'});
	$pepfh->write_seq($pepo);
    }
    $seqfh->write_seq($seqo);
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
