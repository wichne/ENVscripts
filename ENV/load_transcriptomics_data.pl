#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

my %args;
&getopts('D:u:p:f:', \%args);

my $dbh = ENV::connect(\%args);

my $filename = $args{f};
if (! $filename) { die "You need to specify a tsv transcriptomics infile with -f.\n"; }
open (my $file, $filename) or die "Can't open file '$file': $!\n";
my $expt = $filename;
$expt =~ s/\..*//;

my $headerrow = <$file>;
chomp $headerrow;
my $accrow;
my @col;
my @h = split(/\t/, $headerrow);
# Parse the header line for where geneid and counts are
for (my $i=0;$i<@h;$i++) {
    if ($h[$i] =~ /count/i) {
	$col[$i]= $h[$i];
    } elsif ($h[$i] =~ /gene-?id/i) {
	$accrow = $i;
    }
}

my %DATA;
my %TOT;
my %SET;
my $q = "SELECT s.set_id, m.feature_id FROM sequence_sets s, seq_set_link l, seq_feat_mappings m"
    . " WHERE l.seq_id=m.seq_id"
    . " AND s.set_id=l.set_id"
    . " AND s.is_current = 1";
my $r = $dbh->selectall_arrayref($q);
foreach my $row (@$r) {
    my ($set_id, $fid) = @$row;
    $SET{$fid} = $set_id;
}

while (my $row = <$file>) {
    chomp $row;
    my @d = split/\t/, $row;

    # from the accession row, grab the accession and the feature_id
    my $acc = $d[$accrow];
    my $fid;
    my $set_id;
    if ($acc =~ /\|/) {
	my @a = split/\|/, $acc;
	foreach $ac(@a) {
	    if ($ac =~ /gnl\_/) {
		$ac =~ s/gnl\_//;
	    }
	    $fid = get_feature_id_by_accession($dbh, $ac);
	    if ($fid) { last }
	}
    } else {
	if ($acc =~ /gnl\_/) {
	    $acc =~ s/gnl\_//;
	}
	$fid = get_feature_id_by_accession($dbh, $acc);
    }
    if (! $fid) {
	warn "$acc: couldn't find a feature_id for it.\n";
	next;
    }

    # drink in all the info and tally the count/set/sample
    for (my $i = 0; $i<@d; $i++) {
	if ($i == $accrow) { next }
	elsif ($col[$i] && $d[$i]) {
	    $DATA{$fid}->{$expt}->{$col[$i]} = $d[$i];
	    $TOT{$SET{$fid}}->{$expt}->{$col[$i]} += $d[$i];
	}
    }
}

foreach my $fid (sort {$SET{$a} <=> $SET{$b} || $a<=>$b } keys %DATA) {
    foreach my $expt (keys %{$DATA{$fid}}) {
	foreach my $samp (sort keys %{$DATA{$fid}->{$expt}}) {
	    my $rel_count = sprintf "%.2f", ($DATA{$fid}->{$expt}->{$samp}*100/$TOT{$SET{$fid}}->{$expt}->{$samp});
#	    print "$SET{$fid}\t$fid\t$expt\t$samp\t" 
#		. $DATA{$fid}->{$expt}->{$samp} 
#	    . "\t" . $TOT{$SET{$fid}}->{$expt}->{$samp}
#	    . "\t$rel_count\n";
	    &insert_transcriptomics_data($dbh, $fid, $expt, $samp, $DATA{$fid}->{$expt}->{$samp}, $rel_count);
#	    &add_abs_value($dbh, $fid, $expt, $samp, $DATA{$fid}->{$expt}->{$samp});
	}
    }
}

sub insert_transcriptomics_data {
    my $dbh = shift;
    my ($fid, $expt, $samp, $abs_val, $norm_val) = @_;

    my $q = "INSERT INTO transcriptomics_data (feature_id, experiment, sample, abs_value, norm_value) VALUES (?, ?, ?, ?, ?)";
    
    my $ret = $dbh->do($q, {}, ($fid, $expt, $samp, $abs_val, $norm_val));
    if (! defined $ret) { print STDERR "oops on $fid, $expt, $samp, $abs_val, $norm_val: $q\n"; }
}

sub add_abs_value {
    my $dbh = shift;
    my ($fid, $expt, $samp, $abs_val) = @_;

    my $q = "UPDATE transcriptomics_data set abs_value=? where feature_id=? and experiment=? and sample=?";
    
    my $ret = $dbh->do($q, {}, ($abs_val, $fid, $expt, $samp));
    if (! defined $ret) { print STDERR "oops2 on $fid, $expt, $samp, $abs_val : $q\n"; }
}


