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
    open($OUT, ">$outfile") || die "Can't open $outfile for write:$!";
} else {
    $OUT = \*STDOUT;
}

# Keep track of this for IGV
my %GENE_IDX;

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

#my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");
my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS", "tRNA", "rRNA", "ncRNA");

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {

    # The sequence accession
    my $saccref = &get_accession_by_seq_id($dbh, $protref->{$fid}->{'location'}->{'seq_id'});
    my $sacc = $saccref->[0];

    # The location stuff
    my ($min, $max) = ($protref->{$fid}->{'location'}->{'feat_min'},
		       $protref->{$fid}->{'location'}->{'feat_max'});
    my $strand = $protref->{$fid}->{'location'}->{'strand'};
    if ($strand eq "+" || $strand == 1) { $strand = "+" }
    elsif ($strand eq "-" || $strand == -1) { $strand = "-" }
    elsif ($strand eq "0") { $strand = "." }
    else { $strand = "?" }

    my $phase = $protref->{$fid}->{'location'}->{'phase'};
    if (! defined $phase) { warn "No phase reported for CDS $fid: defaulting to phase=0\n";
			    $phase = 0; }
    elsif ($phase !~ /^[012]$/) { warn "Invalid phase ('$phase') for CDS $fid: defaulting to phase=0\n";
				  $phase = 0; }

    # The description stuff
    my $ID = "$db|$fid";
    my @alias;
    foreach my $src (keys %{$protref->{$fid}->{'accessions'}}) {
	my $acc = $protref->{$fid}->{'accessions'}->{$src};
#	$ID = "ID=$acc" if ($src eq 'locus_tag' || !$ID);
	push @alias, "$src=$acc";
    }
    my @desc = ($ID,"bin=$setname",@alias);
    while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	my ($k) = sort {$a<=>$b} keys %$vref;
	my $val = $vref->{$k}->[0]->{value};
	if ($qual eq "gene") {
	    if (defined $GENE_IDX{$val}) { $val .= "-" . ++$GENE_IDX{$val} }
	    else { $GENE_IDX{$val} = 1 if ($val) }
	}
	push @desc, "$qual=$val";
    }

    print $OUT join("\t", $sacc, "PNNL", $protref->{$fid}->{feat_type},
		    $min, $max, ".", $strand, $phase, join(";",@desc));
    print $OUT "\n";

}
