#!/usr/bin/perl

# Edited by WCN on 3/9/2021

use Getopt::Std;
use DBI;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
$| = 1;
# this is being written initially to handle the detail-tsv output with --report-unannotated
# For this output, only loading the rows that have a '*' in the first position
# This has been changed to work off domtbl output from the hmmer package.

my %arg;
&getopts('D:u:p:i:t:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
my $threshfile = $arg{'t'} ? $arg{t} : "/projects/db/hmms/KOfam/ko_list";
if (! $infile) { die "No KOfamscan file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}

my ($THRESH, $ANNOT) = read_thresholds($threshfile);

open my $in, $infile or die "Can't open $infile: $!\n";

my $ev_d = "DELETE FROM feature_evidence"
    . " WHERE feature_id = ?"
    . " AND ev_type = 'KO'";
my $evdstg = $dbh->prepare($ev_d);

my $ev_i = "INSERT INTO feature_evidence (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
    . " VALUES (?, ?, ?, 'KOfamscan', 'KO', ?, ?, ?, ?, ?)";
my $evstg = $dbh->prepare($ev_i);

my $ann_d = "DELETE FROM feature_annotations"
    . " WHERE source = 'KOfamscan'"
    . " AND feature_id = ?";
my $anndstg = $dbh->prepare($ann_d);

my $ann_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, edit_by, date_edit, ann_rank, source, is_current)"
    . " VALUES (?, ?, ?, CURRENT_USER(), NOW(), 3, \"KOfamscan\", 1)";
my $annstg = $dbh->prepare($ann_i);

# screen domtblout file
my %GOOD_HITS;
while (my $line = <$in>) {
    next if ($line =~ /^#/);
    chomp $line;

    # Input is space-delimited, and should be 23 fields domtblout
    my @f = split(/\s+/, $line, 23);
    
    my @acc = split/\|/, $f[0];
    my $fid;
    foreach my $acc (@acc) {
	$fid = &get_feature_id_by_accession($dbh, $acc);
	last if $fid;
    }
    if (!$fid) {die "Why no fid for $f[0]?";}
    
    # hit must meet threshold (if there is a threshold)
    # if there's no threshold if it's better than 1e-10
    # hashes indicate comment lines
    if (defined $THRESH->{$f[3]}) { # check score against threshold
	next if ($f[7] < $THRESH->{$f[3]})
    } elsif (! defined $GOOD_HITS{$fid} && $f[6] > 1e-25) {
	next;
    } else {
	#print ".";
    }
    
    
    #my $len = &get_feature_product_length($dbh, $fid);
    #print "$f[3]\t$f[0]\t$THRESH{$f[3]}\t$f[7]\t$f[6]\n";
    push @{$GOOD_HITS{$fid}->{$f[3]}}, \@f;
}

# load rows
foreach my $fid (keys %GOOD_HITS) {
    # delete existing KO ev and ann for fid
    $evdstg->execute($fid);
    $anndstg->execute($fid);
    # find the best hit
    my $bestScore;
    my $bestKO;
    foreach my $ko (keys %{$GOOD_HITS{$fid}}) {
	# the total score is the same for all domain listings, so just look at the first
	my $thisScore = $GOOD_HITS{$fid}->{$ko}->[0]->[7];
	# print "$ko == $thisScore (/$THRESH->{$ko})\n";
	if ($thisScore > $bestScore) {
	    $bestScore = $thisScore;
	    $bestKO = $ko;
	}
    }

    # load the best KO ev domains
    foreach my $f (@{$GOOD_HITS{$fid}->{$bestKO}}) {
	# insert each good hit
	$evstg->execute($fid, $f->[18], $f->[19], $f->[3], $f->[16], $f->[17], $f->[2], $f->[7]);
    }
    # load annotation row
    $annstg->execute($fid, 66, $ANNOT->{$bestKO}->{product});
    if (defined $ANNOT->{$bestKO}->{EC}) {
	$annstg->execute($fid, 1, $ANNOT->{$bestKO}->{EC});
    }
		   
}

sub get_feature_product_length {
    my $dbh = shift;
    my $fid = shift;
    my $fq = "SELECT length(product) FROM sequence_features WHERE feature_id=$fid";
    my @l = $dbh->selectrow_array($fq);
    if ($l[0]==0) {
	my $mq = "SELECT (feat_max - feat_min + 1)/3 from seq_feat_mappings m where feature_id=$fid";
	@l = $dbh->selectrow_array($mq);
    }
    return $l[0];
}

sub read_thresholds {
    my $threshfile = shift;
    my %THRESH;
    my %ANNOT;
    open(my $threshfh, $threshfile) or die "Can't open threshold file $threshfile: $!\n";
    while (my $line = <$threshfh>) {
	chomp $line;
	my ($knum,
	    $threshold,
	    $score_type,
	    $profile_type,
	    $Fmeasure,
	    $nseq,
	    $nseq_used,
	    $alen,
	    $mlen,
	    $eff_nseq,
	    $re_pos,
	    $definition) = split(/\t/, $line);
	$THRESH{$knum} = $threshold unless ($threshold !~ /^\d+\.\d+$/);
	if ($definition =~ /(.*)( \[EC:(.*)\])?/) {
	    $ANNOT{$knum}->{product} = $1;
	    $ANNOT{$knum}->{EC} = $2;
	} else {
	    $ANNOT{$knum}->{product} = $definition
	}
    }
    close($threshfh);
    return \%THRESH, \%ANNOT;
}
