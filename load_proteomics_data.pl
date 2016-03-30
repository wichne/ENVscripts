#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

my %args;
&getopts('D:u:p:f:d', \%args);

my $dbh = ENV::connect(\%args);
my $DEBUG = $args{'d'};

my $filename = $args{f};
if (! $filename) { die "You need to specify a csv proteomics infile with -f.\n"; }
open (my $file, $filename) || die "Can't open file '$file': $!\n";
my $expt = $filename;
$expt =~ s/\..*//;

my $headerrow = <$file>;
chomp $headerrow;
my $accrow;
my @col;
my @h = split(/\t/, $headerrow);
for (my $i=0;$i<@h;$i++) {
    if ($h[$i] =~ /(.*)[\.\s]*(filtered|normalized) spectra/i) {
	my $sample = $1 ? $1 : 'Sample';
	$col[$i]= $sample;
	print STDERR "Found sample $sample in column $i\n" if $DEBUG;
    } elsif ($h[$i] =~ /accession/i) {
	$accrow = $i;
	print STDERR "Accessions are in column $i\n" if $DEBUG;
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
	    print STDERR "Found feature_id $fid\n" if ($DEBUG);
	} else {
	    my @a = split/\|/, $acc;
	    foreach $ac(@a) {
		$fid = get_feature_id_by_accession($dbh, $ac);
		print STDERR "Found feature_id $fid for $ac\n" if ($DEBUG);
		if ($fid) { last }
	    }
	}
    } else {
	$fid = get_feature_id_by_accession($dbh, $acc);
	print STDERR "Found feature_id $fid for $acc\n" if ($DEBUG);
    }
    if (! $fid) {
	warn "Couldn't load $acc because I couldn't find a feature_id for it.\n";
	next;
    }

    for (my $i = 0; $i<@d; $i++) {
	if ($i == $accrow) { next }
	elsif (defined $col[$i] && defined $d[$i]) {
#	    &insert_proteomics_data($dbh, $fid, $expt, $col[$i], $d[$i]);
	    &update_proteomics_data($dbh, $fid, $expt, $col[$i], $d[$i]);
	}
    }
}

sub insert_proteomics_data {
    my $dbh = shift;
    my ($fid, $expt, $samp, $val) = @_;

#    my $q = "INSERT INTO proteomics_data (feature_id, experiment, sample, filt_spectra) VALUES (?, ?, ?, ?)";
    my $q = "UPDATE proteomics_data set norm_value=? where feature_id=? and experiment=? and sample=?";
    if ($DEBUG) {
	print STDERR "\t$fid, $expt, $samp, $val\n";
    } else {
	my $ret = $dbh->do($q, {}, ($fid, $expt, $samp, $val));
	if (! defined $ret) { print STDERR "oops on $fid, $expt, $samp, $val: $q\n"; }
    }
}
sub update_proteomics_data {
    my $dbh = shift;
    my ($fid, $expt, $samp, $val) = @_;
    if ($val == 0) { return }

#    my $q = "INSERT INTO proteomics_data (feature_id, experiment, sample, filt_spectra) VALUES (?, ?, ?, ?)";
    my $q = "UPDATE proteomics_data set norm_value=? where feature_id=? and experiment=? and sample=?";
    if ($DEBUG) {
	print STDERR "\t$fid, $expt, $samp, $val\n";
    } else {
#	my $ret = $dbh->do($q, {}, ($fid, $expt, $samp, $val));
	my $ret = $dbh->do($q, {}, ($val, $fid, $expt, $samp));
	if (! defined $ret) { print STDERR "oops on $fid, $expt, $samp, $val: $q\n"; }
    }
}
