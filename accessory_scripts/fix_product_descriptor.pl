#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;

my %arg;
&getopts('D:u:p:', \%arg);
my $db = $arg{'D'};
my $dbh = &connect(\%arg);

my $q = "SELECT feature_id, value, ann_rank, source"
    . " FROM feature_annotations"
    . " WHERE data_type_id=66";

my $result = $dbh->selectall_arrayref($q);

my $val_u = "UPDATE feature_annotations set value = ?"
    . " WHERE feature_id=?"
    . " AND value = ?";
my $val_sth = $dbh->prepare($val_u);

foreach my $row(@$result) {
    my ($fid, $val, $rank, $source) = @$row;

    if ($val eq "hypothetical protein") { next }
    
    my @ecs;
    my $new_val = $val;

    if ($new_val =~ /^Similar to (.*)/i) {
	$new_val = "$1-related protein";
    }
    
    if ($new_val =~ /\bprobable\b/i) {
	$new_val =~ s/\bprobable\b/putative/ig;
    }

    if ($new_val =~ /\band related\b\b/) {
	$new_val =~ s/\s*and related.*/ family protein/;
    }
    
    if ($new_val =~ /\b(COG\d{4}\:)/i) {
	$new_val =~ s/$1\s*//i;
    }

    if ($new_val =~ /\b(FIG\d{6,8}\:)/i) {
	$new_val =~ s/$1\s*//i;
    }

    if ($new_val =~ /^(FOG\:)/) {
	$new_val =~ s/^$1\s*//;
    }

    if ($new_val =~ /\b(IPR\d{6}\:)/) {
	$new_val =~ s/$1\s*//;
    }

    if ($new_val =~ /\b(SCO\d{4})/) {
	$new_val =~ s/$1\s*//;
    }

    if ($new_val =~ /\b(UPF\d{4})/) {
	$new_val =~ s/$1\s*//;
    }

    if ($new_val =~ /\;/) {
	$new_val =~ s/\s+\;\s+/\//g;
    }

    if ($new_val =~ /(\s*\(?(possibly|probably).*)/i) {
	$new_val =~ s/\Q$1\E//i;
    }

    if ($new_val =~ /sulphide/i) {
	$new_val =~ s/sulphide/sulfide/ig;
    }

    if ($new_val =~ /isation/) {
	$new_val =~ s/isation/ization/g;
    }

    if ($new_val =~ /haem/i) {
	$new_val =~ s/haem/heme/ig;
    }

    if ($new_val =~ /^domain of unknown function/i) {
	$new_val =~ s/^domain/protein/i;
    }

    if ($new_val =~ /^hypothetical protein .+/i) {
	$new_val = "hypothetical protein";
    }

    if ($new_val =~ /^(hypothetical|probable|possible|predicted|uncharacterized)/i) {
	$new_val =~ s/^$1/putative/i;
    }

    my $misspellings = {'methly' => 'methyl',
			    'reducatse' => 'reductase',
			    'diacyglycerol' => 'diacylglycerol',
			    'hydolase' => 'hydrolase',
			    'hypthetical' => 'hypothetical',
			    'protien' => 'protein',
			    'protei\b' => 'protein',
			    'thiamin\/thiamin\b' => 'thiamin/thiamine',
    };

    while (my ($bad, $good) = each %$misspellings) { 
	if ($new_val =~ /$bad/i) {
	    $new_val =~ s/$bad/$good/ig;
	}
    }

    if ($new_val =~ /protein protein/i) {
	$new_val =~ s/protein protein/protein/i;
    }

    if ($new_val =~ /(\s*homolog)/i) {
	$new_val =~ s/$1/\-like/i;
    }

    if ($new_val =~ /(\s*\(?fragment\)?)/i) {
	$new_val =~ s/$1//i;
    }

    if ($new_val =~ /\b(truncat\w+)/i) {
	$new_val =~ s/\s*$1//i;
    }

    if ($new_val =~ /^miscellaneous/i) {
	$new_val = "hypothetical protein";
    }

    $new_val =~ s/^\s+|\s+$//g;
    
    if ($new_val eq "expressed protein"
	|| $new_val eq "conserved protein"
	|| $new_val eq "conserved hypothetical protein"
	|| $new_val eq "protein"
	|| $new_val eq "protein conserved in bacteria"
	|| $new_val =~ /^similar\b/) {
	$new_val = "hypothetical protein";
    }
    
    if ($new_val ne $val) {
	print "$fid :: '$val' ->\n";
	print "\t'$new_val'\n";
	$val_sth->execute($new_val, $fid, $val);
    }
}
