#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

my $opt = {};
&getopts('D:i:n:o:A:', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};
my $acc_source = $opt->{A};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'mySQL@cce55');


my $OUT;
if ($outfile) {
    open $OUT, ">$outfile" || die "Can't open $outfile for write:$!";
} else {
    $OUT = \*STDOUT;
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }
if (!$setname && $setid) { $setname = &set_id_to_set_name($dbh, $setid); }
my $protref = &get_seq_features_by_set_id($dbh, $setid);

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
    my $ID = "ID=${db}_$fid";
    if ($acc_source) {
	$ID = "ID=" . $protref->{$fid}->{'accessions'}->{$acc_source};
    }
    my @alias;
    foreach my $src (keys %{$protref->{$fid}->{'accessions'}}) {
	my $acc = $protref->{$fid}->{'accessions'}->{$src};
#	$ID = "ID=$acc" if ($src eq 'locus_tag' || !$ID);
	push @alias, "$src=$acc";
    }
    my @desc = ($ID,"bin=$setname",@alias);
    while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	my ($k) = sort {$a<=>$b} keys %$vref;
	push @desc, "$qual=$vref->{$k}->[0]->{value}";
    }

    print $OUT join("\t", $sacc, "PNNL", $protref->{$fid}->{feat_type},
		    $min, $max, ".", $strand, $phase, join(";",@desc));
    print $OUT "\n";

}
