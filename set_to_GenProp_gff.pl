#!/usr/bin/perl

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

my $opt = {};
&getopts('D:i:n:o:', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'access');


my $OUT;
if ($outfile) {
    open $OUT, ">$outfile" || die "Can't open $outfile for write:$!";
} else {
    $OUT = \*STDOUT;
}
print $OUT "##gff-version 3\n";

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {

    # The sequence accession
    my $saccref = &get_accession_by_seq_id($dbh, $protref->{$fid}->{'location'}->{'seq_id'});
    my $sacc = $saccref->[0];

    # The feature type
    my $feat_type = $protref->{$fid}->{feat_type};
    if ($feat_type =~ /^[mnorst]*RNA$/) {
	if ($feat_type ne "tRNA") { $feat_type = "RNA" }
    }

    # The location stuff
    my ($min, $max) = ($protref->{$fid}->{'location'}->{'feat_min'},
		       $protref->{$fid}->{'location'}->{'feat_max'});
    my $strand = $protref->{$fid}->{'location'}->{'strand'};
    if ($strand eq "+" || $strand == 1) { $strand = "+" }
    elsif ($strand eq "-" || $strand == -1) { $strand = "-" }
    elsif ($strand eq "0" || $strand eq ".") { $strand = "." }
    else { $strand = "?" }

    my $phase = $protref->{$fid}->{'location'}->{'phase'};
    if (! defined $phase) { warn "No phase reported for CDS $fid: defaulting to phase=0\n";
			    $phase = 0; }
    elsif ($phase !~ /^[012]$/) { warn "Invalid phase ('$phase') for CDS $fid: defaulting to phase=0\n";
				$phase = 0; }

    # The description stuff
    my @accs;
    foreach my $src (sort keys %{$protref->{$fid}->{'accessions'}}) {
	push @accs, $protref->{$fid}->{'accessions'}->{$src};
    }
    my $acc = shift @accs;
    my $alias = join(",",@accs);

    # if it's a CDS, print a gene feature...
    if ($feat_type eq "CDS") {
	my @desc = ("ID=gene$fid", "locus_tag=$acc", "bin=$setname");
	push @desc, "alias=$alias" if ($alias);
#	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
#	    my ($k) = sort {$a<=>$b} keys %$vref;
#	    push @desc, "$qual=\"$vref->{$k}->{value}\"";
#	}
	
	print $OUT join("\t", $sacc, "PNNL", "gene",
			$min, $max, ".", $strand, ".", join(";",@desc));
	print $OUT "\n";

	# now print the CDS feature
	my @desc = ("ID=$feat_type$fid", "Parent=gene$fid", "feat_name=$acc", "locus_tag=$acc", "bin=$setname");
#       push @desc, "alias=$alias" if ($alias);
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
#	    my ($toprank) = sort {$a<=>$b} keys %$rref;
	    $qual = lc($qual);
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    my $val = $vref->{$k}->[0]->{value};
	    $val =~ s/\,/\%2C/g;
	    $val =~ s/\;/\%3B/g;
	    $val =~ s/\=/\%3D/g;
	    push @desc, "$qual=$val";
	}
    
	print $OUT join("\t", $sacc, "PNNL", $feat_type,
			$min, $max, ".", $strand, $phase, join(";",@desc));
	print $OUT "\n";
    } else {
	my @desc = ("ID=$feat_type$fid", "feat_name=$acc", "locus_tag=$acc", "bin=$setname");
	push @desc, "alias=$alias" if ($alias);
	while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	    $qual = lc($qual);
	    my ($k) = sort {$a<=>$b} keys %$vref;
	    my $val = $vref->{$k}->[0]->{value};
	    $val =~ s/\,/\%2C/g;
	    $val =~ s/\;/\%3B/g;
	    $val =~ s/\=/\%3D/g;
	    push @desc, "$qual=$val";
	}
    	print $OUT join("\t", $sacc, "PNNL", $feat_type,
			$min, $max, ".", $strand, $phase, join(";",@desc));
	print $OUT "\n";
    }

    # Now print any HMM evidence
    my $hmm_ev_ref = &get_evidence($dbh, $fid, "HMM");
    foreach my $hmm_ref (@$hmm_ev_ref) {
	# calculate the genome position
	my ($ev_min, $ev_max);
	if ($strand eq "-") {
	    $ev_min = $max - ($hmm_ref->[2] * 3 - 3);
	    $ev_max = $max - ($hmm_ref->[1] * 3 - 3);
	} else {
	    $ev_min = $min + ($hmm_ref->[1] * 3 - 3);
	    $ev_max = $min + ($hmm_ref->[2] * 3 - 3);
	}
	if ($ev_min < $min || $ev_max > $max) { die "Feature $fid $acc ($min/$max, $strand) has evidence for $hmm_ref->[4] which falls outside gene range ($hmm_ref->[1]/$hmm_ref->[2] -> $ev_min/$ev_max).\n"; }
	my @desc = ("ID=$hmm_ref->[4].$acc", "Parent=CDS$fid",
		    "Target=$hmm_ref->[4] $hmm_ref->[1] $hmm_ref->[2]",
		    "trusted_cutoff=$hmm_ref->[5]", "noise_cutoff=$hmm_ref->[6]");
	print $OUT join("\t", $sacc, "PNNL", "HMM_hit",
			$ev_min, $ev_max, $hmm_ref->[3], $strand, ".", join(";",@desc));
	print $OUT "\n";
    }
}

sub get_evidence {
    my $dbh = shift;
    my $feature_id = shift;
    my $ev_type = shift;
    my $hmm_q = "SELECT e.feature_id, e.feat_min, e.feat_max, e.score, "
	. " h.hmm_acc, h.trusted_cutoff, h.noise_cutoff"
	. " FROM feature_evidence e, egad.hmm2 h"
	. " WHERE e.ev_accession=h.hmm_acc"
	. " AND e.feature_id=$feature_id"
	. " AND e.ev_type='$ev_type'"
	. " ORDER BY e.score DESC";
    my $ref = $dbh->selectall_arrayref($hmm_q);
    return $ref;
}
