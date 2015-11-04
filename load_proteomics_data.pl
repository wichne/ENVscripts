#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

my %args;
&getopts('D:u:p:f:', \%args);

my $dbh = ENV::connect(\%args);

my $filename = $args{f};
if (! $filename) { die "You need to specify a csv proteomics infile with -f.\n"; }
open (my $file, $filename) or die "Can't open file '$file': $!\n";
my $expt = $filename;
$expt =~ s/\..*//;

my $headerrow = <$file>;
chomp $headerrow;
my $accrow;
my @col;
my @h = split(/\t/, $headerrow);
for (my $i=0;$i<@h;$i++) {
    if ($h[$i] =~ /(.*)\.?filtered spectra/i) {
	my $sample = $1 ? $1 : 'Sample';
	$col[$i]= $sample;
    } elsif ($h[$i] =~ /accession/i) {
	$accrow = $i;
    }
}

while (my $row = <$file>) {
    chomp $row;
    my @d = split/\t/, $row;
    my $acc = $d[$accrow];
    my $fid;
    if ($acc =~ /\|/) {
	if ($acc =~ /hotlake_ucc\\|(\d+)$/) {
	    $fid = $1;
	} else {
	    my @a = split/\|/, $acc;
	    foreach $ac(@a) {
		$fid = get_feature_id_by_accession($dbh, $ac);
		if ($fid) { last }
	    }
	}
    } else {
	$fid = get_feature_id_by_accession($dbh, $acc);
    }
    if (! $fid) {
	warn "Couldn't load $acc because I couldn't find a feature_id for it.\n";
	next;
    }

    for (my $i = 0; $i<@d; $i++) {
	if ($i == $accrow) { next }
	elsif ($col[$i] && $d[$i]) {
	    &insert_proteomics_data($dbh, $fid, $expt, $col[$i], $d[$i]);
	}
    }
}

sub insert_proteomics_data {
    my $dbh = shift;
    my ($fid, $expt, $samp, $val) = @_;

    my $q = "INSERT INTO proteomics_data (feature_id, experiment, sample, filt_spectra) VALUES (?, ?, ?, ?)";
    
    my $ret = $dbh->do($q, {}, ($fid, $expt, $samp, $val));
    if (! defined $ret) { print STDERR "oops on $fid, $expt, $samp, $val: $q\n"; }
}
