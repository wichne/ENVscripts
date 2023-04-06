#!/usr/bin/perl

use strict;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use lib $ENV{ENVSCRIPTS};
#use ProkGene;
use ENV;

my %opt;
&getopts('D:i:u:p:t:', \%opt);

if ($opt{h}) {
    die "USAGE\nset_to_pep.pl -D db -u user -p pswd [ -i set_id(all) -t translation_table(11) ]\n";
}

# specify translation table
my $tt = 11;
if ($opt{'t'}) { $tt = $opt{'t'} }

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = &connect(\%opt);

my @setid = $opt{i} ? $opt{i} : ();

if (!@setid) { @setid = &set_ids($dbh) }

foreach my $setid (@setid) {
    print STDERR "Doing set $setid...\n";
	my $seq_ids = set_id_to_seq_ids($dbh, $setid);
	foreach my $seq_id(@$seq_ids) {
		print STDERR "Doing sequence $seq_id\n";
	    my $protref = get_features_by_seq_id($dbh, $seq_id);
	    my $seqref =  get_sequence_by_seq_id($dbh, $seq_id);
	    foreach my $fid (keys %$protref) {
			if ($protref->{$fid}->{'feat_type'} ne "CDS") { next }
			my ($min, $max, $strand,
				$min_partial, $max_partial, $phase) = (
								$protref->{$fid}->{'location'}->{$seq_id}->{'feat_min'},
								$protref->{$fid}->{'location'}->{$seq_id}->{'feat_max'},
								$protref->{$fid}->{'location'}->{$seq_id}->{'strand'},
								$protref->{$fid}->{'location'}->{$seq_id}->{'min_partial'},
								$protref->{$fid}->{'location'}->{$seq_id}->{'max_partial'},
								$protref->{$fid}->{'location'}->{$seq_id}->{'phase'});
			if (!($min && $max) || ($min >= $max)) { warn "Why $min/$max for $fid?\n";}
			my $seq = &get_subseq($seqref, $seq_id, $min, $max, $strand);
			my $seqo = Bio::Seq->new(-display_id => $fid,
						-seq => $seq);
			my $complete = $min_partial || $max_partial ? 0 : 1; 
			my $protobj = $seqo->translate(-complete => $complete,
							-frame    => $phase,
							-codontable_id => $tt);
			&update_product($dbh, $fid, $protobj->seq);
		}
    }
}

sub get_subseq {
    my ($SEQ, $seq_id, $min, $max, $strand) = @_;
#	this would be for a BioPerl Seq obj
    my $subseq = $SEQ->{$seq_id}->trunc($min, $max);
#	but what we have is a hash with sequence data, so:
#	my $subseq = substr($SEQ->{$seq_id}->{sequence}, $min - 1, $max - $min + 1);
    if ($strand < 0 || $strand eq "-") {
		return $subseq->revcom->seq;
    } else {
		return $subseq->seq;
    }
}
