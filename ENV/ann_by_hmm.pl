#!/usr/bin/perl

#08/21/2016 JMo Edited to take in feature_id ($acc[1]) instead of using subroutine get_feature_by_accession for hmm loading

use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use strict;

my %arg;
&getopts('D:i:u:p:s:', \%arg);

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $user = $arg{u} ? $arg{u} : $ENV{USER};
my $password = $arg{p} or die "No password provided.";
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg{D}", $user, $password);

my $hmm_file = $arg{i};

my $hmm_q = "SELECT * FROM egad.hmm2";
my $hmmref = $dbh->selectall_hashref($hmm_q, 'hmm_acc');

my $set_id = $arg{s};

if ($hmm_file) {
    print STDERR "Loading data from $hmm_file...\n";
    open my $in, $hmm_file or die "Can't open $hmm_file: $!\n";
    
    ## 1. Read in line
    while (my $line = <$in>) {
	next if ($line =~ /^#/);
	chomp $line;
	my ($target,
	    $t_acc,
	    $tlen,
	    $q_acc,
	    $q_name,
	    $qlen,
	    $full_evalue,
	    $full_score,
	    $full_bias,
	    $dom_index,
	    $dom_tot,
	    $c_evalue,
	    $i_evalue,
	    $dom_score,
	    $dom_bias,
	    $hmm_lo,
	    $hmm_hi,
	    $ali_lo,
	    $ali_hi,
	    $env_lo,
	    $env_hi,
	    $accuracy,
	    $target_desc) = split/\s+/, $line;

	## 2. Grab the prot acc and the hmm acc

	# all I really care about is the prot_acc and the hmm_acc.
	# the rest comes from the db.
	my $hmm_acc;
	my @acc;
	if ($q_acc =~ /^TIGR\d{5}/ || $q_acc =~ /^PF\d{5}/) {
	    $hmm_acc = $q_acc;
	    @acc = split(/\|/, $target);
	} elsif ($t_acc =~ /^TIGR\d{5}/ || $t_acc =~ /^PF\d{5}/) {
	    $hmm_acc = $t_acc;
	    @acc = split(/\|/, $q_acc);
	} else {
	    die "Confused about which field contains hmm_acc - time for some software engineering!\n";
	}

	&load_hmm_annotation($dbh, $acc[1], $hmm_acc)
	    if ($hmmref->{$hmm_acc}->{iso_type} eq "equivalog" &&
		$full_score > $hmmref->{$hmm_acc}->{'trusted_cutoff'});
    }

} else {
    # Get hmm hit info from database
    my @set_ids;
    if ($arg{s}) { push @set_ids, $arg{s} }
    else { @set_ids = &set_ids($dbh) }
    
    print STDERR "Retrieving HMM hit data from database $arg{D}...\n";
    foreach my $set_id(@set_ids) {
	print STDERR "\tAnnotating set_id : $set_id\n";
	my $feat_count_q = "SELECT count(f.feature_id) FROM sequence_features f"
	    . " INNER JOIN seq_feat_mappings m ON m.feature_id=f.feature_id"
	    . " INNER JOIN seq_set_link l ON m.seq_id = l.seq_id"
	    . " WHERE feat_type='CDS'"
	    . " AND l.set_id = $set_id"
	    ;
	
	my $hmm_feat_count_q = "SELECT count(distinct f.feature_id) FROM feature_evidence f"
	    . " INNER JOIN seq_feat_mappings m ON m.feature_id=f.feature_id"
	    . " INNER JOIN seq_set_link l ON m.seq_id = l.seq_id"
	    . " WHERE ev_type='HMM'"
	    . " AND l.set_id = $set_id"
	    ;
	my ($fc) = $dbh->selectrow_array($feat_count_q);
	my ($hfc) = $dbh->selectrow_array($hmm_feat_count_q);
	print STDERR "I see $hfc features with hmm hits out of $fc features total.\nDoes that sound right?(Y/n) ";
	my $answer = <STDIN>;
	if ($answer =~ /^n/i) { print STDERR " Dying! "; die; }

	my $hmm_q = "SELECT ev.feature_id, ev.ev_accession, ev.score"
	    . " FROM feature_evidence ev"
	    . " INNER JOIN seq_feat_mappings m ON m.feature_id=ev.feature_id"
	    . " INNER JOIN seq_set_link l ON m.seq_id = l.seq_id"
	    . " WHERE ev_type='HMM'"
	    . " AND ev_accession like 'TIGR%'"
	    . " AND l.set_id = $set_id"
	    ;
	#print "$hmm_q\n";
	my $hmm_ev_ref = $dbh->selectall_arrayref($hmm_q);
	foreach my $row (@$hmm_ev_ref) {
	    my $hmm_acc = $row->[1];
	    #print "$hmm_acc\n";
	    &load_hmm_annotation($dbh, $row->[0], $hmm_acc)
		if ($hmmref->{$hmm_acc}->{iso_type} eq "equivalog" &&
		    $row->[2] > $hmmref->{$hmm_acc}->{'trusted_cutoff'});
	}
    }
}

sub load_hmm_annotation {
    my $dbh = shift;
    my $feat_id = shift;
    my $hmm_acc = shift;

    # clear out existing data
    my $del_q = "DELETE FROM feature_annotations"
	. " WHERE feature_id=$feat_id"
	. " AND ann_rank = 3"
	. " AND source = 'TIGRFAM_equivalog'";
    my $del_r = $dbh->do($del_q);
    if (! defined $del_r) {
	die "Error query: $del_q\n";
    }
    
    my $upd_q = "INSERT feature_annotations"
	. " (feature_id, data_type_id, value, ann_rank, source, edit_by)"
	. " VALUES ($feat_id, ?, ?, 3,\"TIGRFAM_equivalog\", USER())"
	;
    my $sth = $dbh->prepare($upd_q);
    if ($hmmref->{$hmm_acc}->{hmm_com_name}) {
	$sth->execute(66, $hmmref->{$hmm_acc}->{hmm_com_name});
    }
    if ($hmmref->{$hmm_acc}->{gene_sym}) {
	$sth->execute(35, $hmmref->{$hmm_acc}->{gene_sym});
	}
    if ($hmmref->{$hmm_acc}->{ec_num}) {
	$sth->execute(1, $hmmref->{$hmm_acc}->{ec_num});
    }
}

exit;
