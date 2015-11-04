#!/usr/bin/perl

use strict;
use Getopt::Std;
use lib $ENV{SCRIPTS};
use ENV;

# output needs to look like 
# start   stop   featurekey
#                              qualifierkey    value

# each CDS will need a gene and a CDS feature
# gene just needs locus_tag

my %args;
&getopt('D:u:p:o:n:s:', \%args);
# D database reqd
# u user def $ENV{USER}
# p password reqd
# o output filename reqd
# n set name reqd
# s source label for accession to be used def 'locus_tag'


my $dbh = connect(\%args);

my $set_name = $args{'n'};
my $set_id = set_name_to_id($dbh, $set_name);
my $source = $args{'s'} ? $args{'s'} : 'locus_tag';

my $output = $args{o} or die "Provide output file with -o\n";
open OUT, ">$output" or die "Can't open '$output' for write: $!\n";

my %GENE_IDX;

my $seqs_ref = get_sequences_by_set_id($dbh, $set_id);
foreach my $sid (sort { $seqs_ref->{$b}->{seq_length} <=> $seqs_ref->{$a}->{seq_length} } keys %$seqs_ref) {
    print OUT ">Feature " . $seqs_ref->{$sid}->{'accessions'}->[0] . "\n";
    my $features_ref = get_features_by_seq_id($dbh, $sid);
    
    my $indent = "\t" x 4;
    select OUT;
	
    foreach my $fid (sort { $features_ref->{$a}->{'location'}->{$set_id}->{'feat_min'} <=> 
				$features_ref->{$b}->{'location'}->{$set_id}->{'feat_min'} ||
				$features_ref->{$a}->{'location'}->{$set_id}->{'feat_max'} <=> 
				$features_ref->{$b}->{'location'}->{$set_id}->{'feat_max'} } keys %$features_ref) {
	my $ref = $features_ref->{$fid};
	if ($sid != $ref->{'location'}->{$set_id}->{'seq_id'}) {
	    #die "Feature $fid: location seq_id (" . $ref->{'location'}->{$set_id}->{'seq_id'} . ") for set $set_id does not match $sid\n";
	    next;
	}
	
#	if ($ref->{'location'}->{$set_id}->{'min_partial'}) {
#	    $ref->{'location'}->{$set_id}->{feat_min} = 
#		$ref->{'location'}->{$set_id}->{'strand'} < 0 ?
#		">" . $ref->{'location'}->{$set_id}->{feat_min} :
#		"<" . $ref->{'location'}->{$set_id}->{feat_min};
#	}
#	if ($ref->{'location'}->{$set_id}->{'max_partial'}) {
#	    $ref->{'location'}->{$set_id}->{feat_max} = 
#		$ref->{'location'}->{$set_id}->{'strand'} < 0 ?
#		"<" . $ref->{'location'}->{$set_id}->{feat_max} :
#		">" . $ref->{'location'}->{$set_id}->{feat_max};
#	}

	my ($end5, $end3);
	if ($ref->{'location'}->{$set_id}->{'strand'} < 0) {
	    $end5 = $ref->{'location'}->{$set_id}->{feat_max};
	    if ($ref->{'location'}->{$set_id}->{'max_partial'}) {
		$end5 = "<" . $end5;
	    }
	    $end3 = $ref->{'location'}->{$set_id}->{feat_min}; 
	    if ($ref->{'location'}->{$set_id}->{'min_partial'}) {
		$end3 = ">" . $end3;
	    }
	} else {
	    $end5 = $ref->{'location'}->{$set_id}->{feat_min};
	    if ($ref->{'location'}->{$set_id}->{'min_partial'}) {
		$end5 = "<" . $end5;
	    }
	    $end3 = $ref->{'location'}->{$set_id}->{feat_max};
	    if ($ref->{'location'}->{$set_id}->{'max_partial'}) {
		$end3 = ">" . $end3;
	    }
	}

	if (!($end5 && $end3)) { die "What is wrong with feature $fid on seq $sid coords $end5/$end3\n";}
	
	my ($br) = sort { $a <=> $b } keys %{$ref->{'annotation'}};
	my ($src) = keys %{$ref->{'annotation'}->{$br}};

	my $gene_sym = defined $ref->{'annotation'}->{$br}->{$src}->{'gene'} ? $ref->{'annotation'}->{$br}->{$src}->{'gene'}->[0] : "";
	if (defined $GENE_IDX{$gene_sym}) { $gene_sym .= "-" . ++$GENE_IDX{$gene_sym} }
	else { $GENE_IDX{$gene_sym} = 1 if ($gene_sym) }

	my $product = "hypothetical protein";
	foreach my $p(@{$ref->{'annotation'}->{$br}->{$src}->{'product'}}) {
	    if ($p ne "") { $product = $p; last; }
	}
	my $feat_type = $ref->{'feat_type'};
	if (grep /$feat_type/, ("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA")) {
	    print "$end5\t$end3\tgene\n";
	    if (!defined $ref->{'accessions'}->{$source} ||
		! @{$ref->{'accessions'}->{$source}}) { die "DIE: Why no $source accession for feature $fid?\n";}
	    else { print $indent, "locus_tag\t" . $ref->{'accessions'}->{$source}->[0] . "\n"; }

	    foreach my $olt (@{$ref->{'accessions'}->{'old_locus_tag'}}) {
		print $indent, "old_locus_tag\t$olt\n";
	    }

	    print $indent, "gene\t" . $gene_sym . "\n" if ($gene_sym);
	    if ($ref->{'location'}->{$set_id}->{'pseudo'} ||
		($ref->{'feat_type'} eq "tRNA" && $product =~ /pseudo/i)) {
		print $indent, "pseudo\n";
		if ($ref->{'feat_type'} eq "tRNA") {
		    $product = "Xxx";
		} else {
		    print $indent, "note\t$product\n";
		}
	    }
	}
	
	if ($ref->{'feat_type'} eq "misc_feature" ||
	    $ref->{'feat_type'} eq "repeat_region") {
	    print "$end5\t$end3\t" . $ref->{'feat_type'}. "\n";
	    print $indent, "note\t$product\n";
	    next;
	}

	if ($ref->{'feat_type'} eq "misc_binding") {
	    print "$end5\t$end3\tmisc_binding\n";
	    print $indent, "note\t$product\n";
	    my $bound_moiety=&get_bound_moiety($dbh, $fid, $product);
	    print $indent, "bound_moiety\t$bound_moiety\n";
	    next;
	}
	
	# do not print out CDS features that are pseudogenes
	if ($ref->{'location'}->{$set_id}->{'pseudo'}) { next }
	# handle ribosomal slippage and selenocysteine CDS
#	elsif (defined $qual_r->{ribosomal_slippage} ||
#	       defined $qual_r->{transl_except} ) {
#	    if (defined $qual_r->{CDS_location}) {
#		my $rev = $qual_r->{CDS_location}->{value} =~ /complement/ ? 1 : 0;
#		my @seg = $rev ?
#		    reverse(split(",",$qual_r->{CDS_location}->{value})) :
#		    split(",",$qual_r->{CDS_location}->{value});
#		for (my $i=0; $i<@seg; $i++) {
#		    $seg[$i] =~ /(\d+)\.{2}(\d+)/;
#		    ($end5, $end3) = $rev ? ($2, $1) : ($1, $2);
#		    print "$end5\t$end3";
#		    if ($i == 0) { print "\t$ref->{feat_type}\n" }
#		    else { print "\n" }
#		}
#	    }
#	} else {
	print "$end5\t$end3\t$ref->{feat_type}\n";
#	}
	
	# this should have the locus_tag in it
	my $locus_tag = $ref->{'accessions'}->{$source}->[0];
	print $indent, "protein_id\tgnl|PNNL|$locus_tag\n";
	if ($ref->{'location'}->{$set_id}->{'phase'}) {
	    my $phase = $ref->{'location'}->{$set_id}->{'phase'} + 1;
	    print $indent, "codon_start\t$phase\n";
	}

	print $indent, "product\t$product\n";
	print $indent, "gene\t" . $gene_sym . "\n" if ($gene_sym);
	if (defined $ref->{'annotation'}->{$br}->{$src}->{'EC_number'}) {
	    foreach my $ec (@{$ref->{'annotation'}->{$br}->{$src}->{'EC_number'}}) {
		my @ec = split(/[\s\,\;]+/, $ec);
		foreach my $e (@ec) {
		    print $indent, "EC_number\t$e\n";
		}
	    }
#	foreach my $q (keys(%$qual_r)) {
#	    if (
#		$q eq "CDS_location" ||
#		$q eq "SO_term"
#		) { next }
#	    elsif ($q eq "EC_number" && $qual_r->{$q}->{value} =~ /\;\;/) {
#		$qual_r->{$q}->{value} =~ s/\"//g;
#		my @e = split /\;\;/, $qual_r->{$q}->{value};
#		foreach my $e (@e) {
#		    print "${indent}$q\t$e\n";
#		}
#	    } else {
#		print "${indent}$q\t$qual_r->{$q}->{value}\n";
#	    }
	}
	if ($feat_type eq "ncRNA") {
	    my $class = &get_ncRNA_class($dbh, $fid, $product);
	    print $indent, "ncRNA_class\t$class\n";
	}
	
    }
}

sub get_bound_moiety {
    my $dbh=shift;
    my $fid = shift;
    my $product = shift;
    if ($product =~ /\b(RF\d{5})\b/) {
	my $rfacc = $1;
	my $q = "SELECT bound_moiety from rfam.bound_moiety where accession=\"$rfacc\"";
	my $r = $dbh->selectcol_arrayref($q);
	return $r->[0];
    } else {
	my $q = "SELECT bound_moiety from rfam.rfam r, rfam.bound_moiety b, feature_evidence e"
	    . " WHERE e.feature_id=$fid AND e.ev_type='RFAM' AND e.ev_accession=b.accession"
	    ;
	my $r = $dbh->selectcol_arrayref($q);
	return $r->[0];
    }
}

sub get_ncRNA_class {
    my $dbh=shift;
    my $fid = shift;
    my $product = shift;
    if ($product =~ /\b(RF\d{5})\b/) {
	my $rfacc = $1;
	my $q = "SELECT ncRNA_class from rfam.rfam where accession=\"$rfacc\"";
	my $r = $dbh->selectcol_arrayref($q);
	return $r->[0];
    } else {
	my $q = "SELECT ncRNA_class from rfam.rfam r, feature_evidence e"
	    . " WHERE e.feature_id = $fid AND e.ev_type='RFAM' AND e.ev_accession=r.accession"
	    ;
	my $r = $dbh->selectcol_arrayref($q);
	return $r->[0];
    }
}
