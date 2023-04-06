#!/usr/bin/perl

#09/16/2014, Jennifer Mobberley (mobberley.jennifer@pnnl.gov)
#This script will load the best-hit dbCAN hmm generated from local annotation scripts (see README file at csbl.bmb.uga.edu/dbCAN/download). An additional step to append dbCAN hmm model length was carried out using the /home/mobb021/hmm/dbCAN/dbCAN_hmm_prep.sh
#02/20/2015 WCN adapted from load_dbCan.pl This will now take a set name as an argument, build a fresh pep file, run the hmmscan (actually hmmsearch) and then load the database with the data (unless you ask it not to)
#08/21/2016 JMo updated dbCAN_hmm_file path and reworked mysql excute to load into feature_evidence based on feature_id which is contained in the sequence headers of the peptide files


use Getopt::Std;
use DBI;
use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;

our $DEBUG = 1;
our $PARSER = "hmmscan-parser.sh";

my %arg;
&getopts('D:u:p:i:n:m:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $set_name = $arg{'n'};
our $dbCAN_hmm_file = $arg{'m'} ? $arg{'m'} : "/projects/db/hmms/dbCAN/V9/dbCAN-HMMdb-V9.txt";



my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
my $dbcan_data;
if (! $infile) {
    if ($set_name) {
	my $build_pep = 1;
	my $pepfile = $set_name . ".pep";
	if (-e $pepfile) {
	    print "Do you want to rebuild $pepfile?(y/N) ";
	    my $a = <STDIN>;
	    unless ($a =~ /^y/i) { $build_pep = 0 }
	}
	$pepfile = &build_pep_file($db, $set_name) if ($build_pep);

	my $search_dbCAN = 1;
	my $hmmfile = $pepfile . ".dbCAN.domtblout";
	if (-e $hmmfile) {
	    print "Do you want to recreate $hmmfile?(y/N) ";
	    my $a = <STDIN>;
	    unless ($a =~ /^y/i) { $search_dbCAN = 0 }
	}
	$hmmfile = &run_dbCAN_hmmsearch($pepfile) if ($search_dbCAN);
	$dbcan_data = &parse_hmmsearch_output($hmmfile);
	print STDERR "Ready to load...\n";
    } else { die "No input file provided with -i or -n\n";}
} else {
    if (! -r $infile) { die "$infile is not readable: $!\n";}
    $dbcan_data = &parse_hmmsearch_output($infile);
    # open my $in, $infile or die "Can't open $infile: $!\n";
    # while (my $line = <$in>) {
    # 	chomp $line;
    # 	my @f = split/\s+/, $line;
    # 	if (@f != 9) { die "Input file $infile is not in dbcan summary format: ". @f . "\n"; }
    # 	my $acc = shift @f;
    # 	push @{$dbcan_data->{$acc}}, \@f;
    # }
}

my $ev_i = "INSERT INTO feature_evidence"
    . " (feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length,score)"
    . " SELECT feature_id , ?, ?, 'dbCAN', 'CAZy', ?, ?, ?, ?, ?"
    . " FROM feature_evidence WHERE feature_id=?";

my $sth = $dbh->prepare($ev_i);
foreach my $x (keys %$dbcan_data) {
    my @acc = split/\|/, $x;
    my $acc = @acc > 1 ? $acc[1] : $acc[0];
    if (! $acc) { die "Why no accession from $x ?\n"; }
    my $ev_d = "DELETE FROM e"
	. " USING feature_evidence e"
	. " WHERE e.feature_id = \"$acc\""
	. " AND e.ev_type='CAZy'";
    $dbh->do($ev_d) unless ($DEBUG);
    foreach my $hitref (@{$dbcan_data->{$x}}) {
	my ($ev_accession,
	    $evalue,
	    $ev_length,
	    $ev_min,
	    $ev_max,
	    $feat_length,
	    $feat_min,
	    $feat_max,
	    $hmm_cover)= @$hitref;
	my $score = "evalue=$evalue;coverage=$hmm_cover;";
	print STDERR "Inserting $acc :: $ev_accession\n" if ($DEBUG);
	$sth->execute($feat_min, $feat_max, $ev_accession, $ev_min, $ev_max, $ev_length, $score, $acc);# unless ($DEBUG);
    }
}

sub build_pep_file {
    my $db = shift;
    my $set_name = shift;
    system("set_to_pep.pl -D $db -n $set_name -o $set_name.pep");
    return "$set_name.pep";
}

sub run_dbCAN_hmmsearch {
    my $pepfile = shift;
    system("hmmsearch --cpu 10 --domtblout $pepfile.dbCAN.domtblout --noali --acc -o $pepfile.dbCAN.out $dbCAN_hmm_file $pepfile");
    return "$pepfile.dbCAN.domtblout";
}

sub parse_hmmsearch_output {
    my $domtbl_out = shift;
    my %dbcanhits;

    # parse the results of the hmmsearch and pipe into the process
    open(my $parseh, "$PARSER $domtbl_out |") or die "Can't use $PARSER on $domtbl_out: $!\n";
    while (my $line = <$parseh>) {
	chomp $line;
	my ($ev_acc,
	    $ev_len,
	    $prot_acc,
	    $prot_len,
	    $evalue,
	    $hmmfrom,
	    $hmmto,
	    $alifrom,
	    $alito,
	    $cov) = split(/\t/,$line);

	# good hit threshold
	if ($evalue < 1e-15 && $cov > 0.35) {
	    push @{$dbcanhits{$prot_acc}}, [$ev_acc,
					    $evalue,
					    $ev_len,
					    $hmmfrom,
					    $hmmto,
					    $prot_len,
					    $alifrom,
					    $alito,
					    $cov];
	}
    }

    # Only take the best assignment (by evalue)
    foreach my $prot_acc(sort keys %dbcanhits) {
	my @hits = sort { $a->[1] <=> $b->[1] } @{$dbcanhits{$prot_acc}};
	my @good;
	foreach my $hit (@hits) {
	    my $good = 1;
	    foreach my $g (@good) {
		# if there's an overlap that is more than 50% of the length of this hit region, it' no good.
		if (&overlaps($hit->[6], $hit->[7], $g->[6], $g->[7])/($hit->[7]-$hit->[6]+1) > 0.5) {
		    $good = 0;
		}
	    }
	    push @good, $hit if ($good);
	}
	$dbcanhits{$prot_acc} = \@good;
    }
    return \%dbcanhits;
}

sub overlaps {
    my @coord = @_;
    my @sort = sort {$a<=>$b} @coord;
    return ((abs($coord[0]-$coord[1])+1)+(abs($coord[2]-$coord[3])+1)-($sort[3]-$sort[0]+1));
}
