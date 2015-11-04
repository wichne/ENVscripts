#!/usr/bin/perl

use lib "/share/scripts/";
use ENV;
use Getopt::Std;
use strict;

my %arg;
&getopts('D:u:p:o:i:n:s:h', \%arg);
my $dbh = ENV::connect(\%arg);

my $usage = "overlap_analysis.pl -D db -p dbpswd [-u dbuser -o outputfile -i set_id -n set_name -s seq_id ]\n";

$| = 1;
if ($arg{h}) { print $usage; exit()}
if ($arg{o}) {
    open my $out, ">$arg{o}";
    select $out;
}

my @set_ids;
if ($arg{s}) {
    @set_ids = split/\,/, $arg{s};
} elsif ($arg{n}) {
    my @set_names = split/\,/, $arg{n};
    foreach my $name (@set_names) {
	push @set_ids, &set_name_to_id($dbh, $name);
    }
} else {
    @set_ids = sort {$a<=>$b} keys(%{&get_seq_sets($dbh)});
}

foreach my $set_id (@set_ids) {
    print STDERR "Checking set $set_id\n";
    my $seq_ids = &set_id_to_seq_ids($dbh, $set_id);
    foreach my $seq_id (@$seq_ids) {
	my $feat_ref = &get_features_by_seq_id($dbh, $seq_id);
	my @feat_ids = sort {$feat_ref->{$a}->{'location'}->{$seq_id}->{'feat_min'} <=>
				 $feat_ref->{$b}->{'location'}->{$seq_id}->{'feat_min'} }
	keys %$feat_ref;
	print STDERR "\tChecking " . scalar(@feat_ids) . " on seq $seq_id\n";
	for (my $i=0; $i<@feat_ids; $i++) {
	    my $sfeat_id = $feat_ids[$i];
	    my $sev_ref = &get_evidence_for_feature($dbh, $sfeat_id);
	    my $slen = $feat_ref->{$sfeat_id}->{'location'}->{$seq_id}->{'feat_max'} - $feat_ref->{$sfeat_id}->{'location'}->{$seq_id}->{'feat_min'} + 1;
	    for (my $j= $i > 10 ? $i - 10 : 0; $j<@feat_ids; $j++) {
		my $ofeat_id = $feat_ids[$j];
		if ($sfeat_id == $ofeat_id) { last }
		my $overlap = &overlap($feat_ref->{$sfeat_id}->{'location'}->{$seq_id}->{'feat_min'},
				       $feat_ref->{$sfeat_id}->{'location'}->{$seq_id}->{'feat_max'},
				       $feat_ref->{$ofeat_id}->{'location'}->{$seq_id}->{'feat_min'},
				       $feat_ref->{$ofeat_id}->{'location'}->{$seq_id}->{'feat_max'});
		my $olen = $feat_ref->{$ofeat_id}->{'location'}->{$seq_id}->{'feat_max'} - $feat_ref->{$ofeat_id}->{'location'}->{$seq_id}->{'feat_min'} + 1;
		if (!$overlap || ($overlap <= 60 && $slen > 120 && $olen > 120)) { next }
		print "$set_id $seq_id $feat_ref->{$sfeat_id}->{feat_type} $sfeat_id ($slen)"
		    . " overlaps $feat_ref->{$ofeat_id}->{feat_type} $ofeat_id ($olen) by $overlap nt. ";
		my $oev_ref = &get_evidence_for_feature($dbh, $ofeat_id);
		if (@$sev_ref == 0 && @$oev_ref > 0 &&
		    $overlap/$slen > 0.5 &&
		    $feat_ref->{$sfeat_id}->{feat_type} eq $feat_ref->{$ofeat_id}->{feat_type}) {
		    print "\t$sfeat_id has no evidence -> hiding $sfeat_id.\n";
		    &not_current($dbh, $sfeat_id);
		} elsif (@$oev_ref == 0 && @$sev_ref > 0 &&
			 $overlap/$olen > 0.5 &&
			 $feat_ref->{$sfeat_id}->{feat_type} eq $feat_ref->{$ofeat_id}->{feat_type}) {
		    print "\t$ofeat_id has no evidence -> hiding $ofeat_id.\n";
		    &not_current($dbh, $ofeat_id);
		} elsif (@$oev_ref > 0 && @$sev_ref > 0 && 
			 ($overlap/$slen > 0.5 || $overlap/$olen > 0.5) ) {
		    print "\tBoth have evidence. Manual resolution required.\n";
		} elsif ($feat_ref->{$sfeat_id}->{feat_type} eq "tRNA" && $overlap/$slen > 0.5 && @$oev_ref == 0 ) {
		    print "\t$ofeat_id has no evidence -> hiding $ofeat_id.\n";
		    &not_current($dbh, $ofeat_id);		    
		} elsif ($feat_ref->{$ofeat_id}->{feat_type} eq "tRNA" && $overlap/$olen > 0.5 && @$sev_ref == 0 ) {
		    print "\t$sfeat_id has no evidence -> hiding $sfeat_id.\n";
		    &not_current($dbh, $sfeat_id);		    
		} else {
		    print "\tManual resolution suggested.\n";
		}
	    }
	}
    }
}

sub overlap {
    my ($x1, $y1, $x2, $y2) = @_;
    my @p = sort {$a<=>$b} ($x1, $y1, $x2, $y2);
    my $o = (abs($y1 - $x1) + 1 + abs($y2-$x2) + 1) - ($p[3] - $p[0] + 1);
    if ($o > 0) { return $o } 
    else { return 0 }
}

sub not_current {
    my $dbh = shift;
    my $feat_id = shift;
    my $q = "UPDATE sequence_features"
	. " SET is_current = 0"
	. " WHERE feature_id=$feat_id";
    $dbh->do($q);
}
