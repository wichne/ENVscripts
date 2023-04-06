#!/usr/bin/perl
#edited by JMo January 12 2016 to update for format of new dbCAN output file after dbCAN parsing
#see readme file in scripts/db/hmm/dbCAN for specifics
use lib $ENV{ENVSCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use DBI;
use Cwd;

my $host = $ENV{'DBSERVER'} ? $ENV{'DBSERVER'} : 'localhost';

my %arg;
&getopts('D:p:u:d:h',\%arg);
my $USER = $arg{'u'} ? $arg{'u'} : $ENV{'USER'}; 
my $PSWD = $arg{'p'} or die "Need to provide database password with -P\n";
my $DB = $arg{'D'} or die "Need to provide database with -D\n";

if ($arg{'h'}) {
    print "USAGE load_dbCAN_to_ENV.pl -D db -P dbpassword -d dbCANoutputdirectorypath [-U user ]\n";
}

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$DB", $USER, $PSWD);

my $path = $arg{'d'} or die "Need to provide dbCAN output directory with -d\n";
cwd($path) or die "Can't change path to $path? : $!\n";

opendir(DIR, $path) or die "cannot open directory: $!\n";
my @files = grep(/\.out$/,readdir(DIR));
foreach my $file (@files) {
    $file = $path . "/" . $file;
    ################################################################################
    ########## dbCAN HMMER
    ################################################################################
    if ($file =~ /hmmer/) {
	open (my $HMM, "$file") or die "Can't open $file. Is path correct? $!\n";

	my $delete_q = "DELETE e FROM feature_evidence e"
	    . " WHERE e.feature_id=?"
	    . " AND e.ev_accession = ?"
	    . " AND e.ev_type = 'CAZy'"
	    . " AND e.program='dbCAN'";
	my $dsth = $dbh->prepare($delete_q);

	my $insert_q = "INSERT INTO feature_evidence"
	    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
	    . " VALUES(?, ?, ?, \"dbCAN\", \"CAZy\", ?, ?, ?, ?, ?)";

	my $isth = $dbh->prepare($insert_q);

	my $linecount = 0;
	while (my $line = <$HMM>) {
	    if ($line =~ /^HMM/) { next }
	    $linecount++;
	    chomp $line;
	    my ($hacc,
		$hmm_len,
		$feat_acc,
		$gene_len,
		$evalue,
		$hmm_start,
		$hmm_end,
		$gene_start,
		$gene_end,
		$cov) = split(/\t/,$line);

	    if ($evalue > 1e-15 || $cov < 0.35) { next }
	    
	    my $fid = &get_feature_id_by_accession($dbh, $feat_acc);
	    if (! $fid) { warn "HMMer: Can't get feature_id from '$line'. Not loading this one.\n"; next } 

	    my $sxs = $dsth->execute($fid, $hacc);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	    my $sxs = $isth->execute($fid, $gene_start, $gene_end, $hacc, $hmm_start, $hmm_end, $hmm_len, $evalue);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	}
    } elsif ($file =~ /HotPep/) {
	################################################################################
	########## HOTPEP/PPR
	################################################################################

	open (my $HPEP, $file) or die "Can't open $file. Is path correct? $!\n";

	my $HP_delete_q = "DELETE e FROM feature_evidence e"
	    . " WHERE e.feature_id = ?"
	    . " AND e.ev_accession = ?"
	    . " AND e.program='HotPep'"
	    . " AND e.ev_type='CAZy'";
	my $HPdsth = $dbh->prepare($HP_delete_q);

	my $HP_insert_q = "INSERT INTO feature_evidence"
	    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, score)"
	    . " VALUES(?, 1, 2, \"HotPep\", \"CAZy\", ?, ?)";

	my $HPisth = $dbh->prepare($HP_insert_q);
	my $linecount;
	while (my $line = <$HPEP>) {
	    if ($line =~ /^CAZy/) { next }
	    $linecount++;
	    chomp $line;
	    my ($cazy,
		$ppr,
		$feat_acc,
		$freq,
		$hits,
		$peps,
		$ECs) = split(/\t/, $line);
	    
	    my $fid = &get_feature_id_by_accession($dbh, $feat_acc);
	    if (! $fid) { warn "HOTPEPPPR: Can't get feature_id from '$line'. Not loading this one.\n"; next } 

	    my $evacc = $cazy . "(" . $ppr . ")";
	    my $sxs = $HPdsth->execute($fid, $evacc);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	    my $sxs = $HPisth->execute($fid, $evacc, $freq);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	}
    } elsif ($file =~ /diamond/) {
	################################################################################
	########## DIAMOND
	################################################################################

	open (my $DIA, $file) or die "Can't open $file. Is path correct? $!\n";

	my $D_delete_q = "DELETE e FROM feature_evidence e"
	    . " WHERE feature_id = ?"
	    . " AND e.ev_accession = ?"
	    . " AND ev_type = 'CAZy'"
	    . " AND program = 'diamond'";
	my $Ddsth = $dbh->prepare($D_delete_q);

	my $D_insert_q = "INSERT INTO feature_evidence"
	    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score)"
	    . " VALUES(?, ?, ?, \"diamond\", \"CAZy\", ?, ?, ?, ?, ?)";

	my $Disth = $dbh->prepare($D_insert_q);
	my $linecount;
	while (my $line = <$DIA>) {
	    if ($line =~ /^Gene/) { next }
	    $linecount++;
	    chomp $line;
	    my ($feat_acc,
		$dcazy,
		$per_id,
		$length,
		$mismatches,
		$gap,
		$gene_start,
		$gene_end,
		$cazy_start,
		$cazy_end,
		$evalue,
		$bits) = split(/\t/, $line);
	    
	    my $fid = &get_feature_id_by_accession($dbh, $feat_acc);
	    if (! $fid) { warn "DIAMOND: Can't get feature_id from '$line'. Not loading this one.\n"; next } 

	    my $cazy;
	    if ($dcazy =~ /((GH|GT|CE|CBM|AA|PL)\d+)/) {
		$cazy = $1;
	    } else { warn "$feat_acc: Can't find CAZy acc in '$dcazy'. Need to update REGEX?\n"; next; }
	    
	    my $sxs = $Ddsth->execute($fid, $cazy);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	    my $sxs = $Disth->execute($fid, $gene_start, $gene_end, $cazy, $cazy_start, $cazy_end,
				      $length, $evalue);
	    if (! $sxs) { print STDERR $dbh->errstr() }
	}
    }
}
