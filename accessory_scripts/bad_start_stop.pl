#!/usr/bin/perl
use Getopt::Std;
use lib $ENV{SCRIPTS};
use ENV;
use strict;

my %arg;
&getopts('D:u:p:f:t', \%arg);
if (! $arg{'D'}) { die "bad_start_stop.pl -D db [-u user] -p pswd -f file [-t (nodbwrite)]\n";}
my $dbh = &connect(\%arg);
my $FILE = Bio::SeqIO->new(-file => $arg{'f'});

my $seqq = "select seq_id, seq_length from sequences where iscurrent=1";
my $SEQ = $dbh->selectall_hashref($seqq, 'seq_id');
my $TEST = $arg{'t'};

while (my $seqr = $FILE->next_seq) {
    my $id = $seqr->primary_id;
    my $seq = $seqr->seq;
    my ($bad_start, $bad_stop, $gap);
    my ($fid, $feat_obj, $locr);

    if ($seq !~ /^[atg]tg/i) {
	$bad_start = 1;
    }
    if ($seq !~ /(taa|tag|tga)$/i) {
	$bad_stop = 1;
    }
    if ($seq =~ /n/i &&
	length($seq) % 3 != 0) {
	$gap = 1;
    }

    if ($bad_start || $bad_stop || $gap) {
	if ($id =~ /\|(\d+)\|/) {
	    $fid = $1;
	    $feat_obj = get_features_by_feature_id($dbh, $fid);
	    if ($feat_obj->{$fid}->{feat_type} ne "CDS") { next }
	    $locr = $feat_obj->{$fid}->{'location'};
	} else { die "Couldn't parse feature id from " . $id . "\n";}
    } else { next }

    if ($gap) {
	print "$id appears to be frame-shifted due to scaffold or sequence gap\n";
	foreach my $set_id (keys %$locr) {
	    my $strand = $locr->{$set_id}->{'strand'} == 1 ? -1 : 1;
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    my $q = "UPDATE seq_feat_mappings SET pseudo=1 where feature_id=$fid and seq_id=$seq_id";
	    if ($TEST) {
		print $q, "\n";
	    } else {
		my $r = $dbh->do($q);
		if (! defined $r) { warn $dbh->errstr }
		next;
	    }
	}
    }

    if ($seq =~ /^(cta|tta|tca).*(ca[cat])$/i) {
	print $id, " appears to be on the wrong strand: ", substr($seq, 0, 3), "..", substr($seq, -3, 3), "\n\n";
	foreach my $set_id (keys %$locr) {
	    my $strand = $locr->{$set_id}->{'strand'} == 1 ? -1 : 1;
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    my $q = "UPDATE seq_feat_mappings SET strand='$strand' where feature_id=$fid and seq_id=$seq_id";
	    if ($TEST) {
		print $q, "\n";
	    } else {
		my $r = $dbh->do($q);
		if (! defined $r) { warn $dbh->errstr }
		next;
	    }
	}
    } elsif ($seq =~ /^(cta|tta|tca)/i) {
	foreach my $set_id (keys %$locr) {
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    if ($locr->{$set_id}->{'feat_max'} >= $SEQ->{$seq_id}->{'seq_length'} - 2 &&
		$locr->{$set_id}->{'strand'} > 0) { 
		print "$id may be on wrong strand: ", substr($seq, 0, 3), "..", substr($seq, -3, 3), "\n\n";
		my $strand = $locr->{$set_id}->{'strand'} == 1 ? -1 : 1;
		my $q = "UPDATE seq_feat_mappings SET strand='$strand', max_partial='1'"
		    . " where feature_id=$fid and seq_id=$seq_id";
		$bad_stop = 0;
		if ($TEST) {
		    print $q, "\n";
		} else {
		    my $r = $dbh->do($q);
		    if (! defined $r) { warn $dbh->errstr }
		    next;
		}
	    }
	}
    } elsif ($seq =~ /(ca[cat])$/i) {
	foreach my $set_id (keys %$locr) {
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    if ($locr->{$set_id}->{'feat_min'} <= 3 &&
		$locr->{$set_id}->{'strand'} > 0) { 
		print "$id may be on wrong strand: ", substr($seq, 0, 3), "..", substr($seq, -3, 3), "\n\n";
		my $strand = $locr->{$set_id}->{'strand'} == 1 ? -1 : 1;
		my $q = "UPDATE seq_feat_mappings SET strand='$strand', min_partial='1'"
		    . " where feature_id=$fid and seq_id=$seq_id";
		$bad_start = 0;
		if ($TEST) {
		    print $q, "\n";
		} else {
		    my $r = $dbh->do($q);
		    if (! defined $r) { warn $dbh->errstr }
		    next;
		}
	    }
	}
    }
    
    if ($bad_start) {
	print $id, " has bad start: ", substr($seq, 0, 3), "\n";
	foreach my $set_id (keys %$locr) {
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    print "$id\tseq:$seq_id ($SEQ->{$seq_id}->{seq_length} nt)\t" . $locr->{$set_id}->{feat_min}
	    . " / " . $locr->{$set_id}->{feat_max} 
	    . "\tstrand:" . $locr->{$set_id}->{strand}
	    . "\n";
	    if ($locr->{$set_id}->{'feat_min'} <= 3 &&
		$locr->{$set_id}->{'strand'} > 0) { 
		my $phase = $locr->{$set_id}->{'feat_max'} % 3;
		print "START:update min_partial for $fid to T, phase to $phase\n\n";
		my $q1 = "update seq_feat_mappings"
		    . " set feat_min = 1, min_partial=1, phase='$phase'"
		    . " where seq_id=$seq_id and feature_id=$fid";
		if ($TEST) { print "$q1\n";}
		else { my $r = $dbh->do($q1);
		       if (! defined $r) { warn $dbh->errstr }}
	    } elsif ($locr->{$set_id}->{'feat_max'} >= $SEQ->{$seq_id}->{'seq_length'} - 2 &&
		$locr->{$set_id}->{'strand'} < 0) { 
		my $phase = (abs($locr->{$set_id}->{'feat_min'} - $locr->{$set_id}->{'feat_max'}) + 1) % 3;
		print "START:update max_partial for $fid to T, phase to $phase\n\n";
		my $q2 = "update seq_feat_mappings"
		    . " set feat_max = " . $SEQ->{$seq_id}->{'seq_length'} . ", max_partial=1, phase = '$phase'"
		    . " where seq_id=$seq_id and feature_id=$fid";
		if ($TEST) { print "$q2\n" }
		else { my $r = $dbh->do($q2);
		       if (! defined $r) { warn $dbh->errstr } }
	    }
	}
    }
    if ($bad_stop) {
	if (! $feat_obj) {
	    if ($id =~ /\|(\d+)\|/) {
		$fid = $1;
		$feat_obj = get_features_by_feature_id($dbh, $fid);
	    } else { die "Couldn't parse feature id from " . $id . "\n";}
	}
	if ($feat_obj->{$fid}->{feat_type} ne "CDS") { next }
	my $locr = $feat_obj->{$fid}->{'location'};
	print $id, " has bad stop: ", substr($seq, -3, 3), "\n";
	foreach my $set_id (keys %$locr) {
	    my $seq_id = $locr->{$set_id}->{'seq_id'};
	    print "$id\t$seq_id($SEQ->{$seq_id}->{seq_length})\t" . $locr->{$set_id}->{feat_min}
	    . "\t" . $locr->{$set_id}->{feat_max} 
	    . "\t" . $locr->{$set_id}->{strand}
	    . "\n";
	    if ($locr->{$set_id}->{'feat_min'} <= 3 &&
		$locr->{$set_id}->{'strand'} < 0) { 
		print "STOP:update min_partial for $fid to T\n\n";
		my $q1 = "update seq_feat_mappings"
		    . " set feat_min=1, min_partial=1"
		    . " where seq_id=$seq_id and feature_id=$fid";
		my $r = $dbh->do($q1) unless $TEST;
		if (!$TEST && ! defined $r) { warn $dbh->errstr }
	    } elsif ($locr->{$set_id}->{'feat_max'} >= $SEQ->{$seq_id}->{'seq_length'} - 2 &&
		$locr->{$set_id}->{'strand'} > 0) { 
		print "STOP:update max_partial for $fid to T\n\n";
		my $q2 = "update seq_feat_mappings"
		    . " set feat_max = " . $SEQ->{$seq_id}->{'seq_length'} . ", max_partial=1"
		    . " where seq_id=$seq_id and feature_id=$fid";
		if ($TEST) { print "$q2\n" } else {
		    my $r = $dbh->do($q2);
		    if (! defined $r) { warn $dbh->errstr } }
	    }
	}
	print "\n";
    }
}

