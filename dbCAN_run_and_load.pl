#!/usr/bin/perl

#09/16/2014, Jennifer Mobberley (mobberley.jennifer@pnnl.gov)
#This script will load the best-hit dbCAN hmm generated from local annotation scripts (see README file at csbl.bmb.uga.edu/dbCAN/download). An additional step to append dbCAN hmm model length was carried out using the /home/mobb021/hmm/dbCAN/dbCAN_hmm_prep.sh
#02/20/2015 WCN adapted from load_dbCan.pl This will now take a set name as an argument, build a fresh pep file, run the hmmscan (actually hmmsearch) and then load the database with the data (unless you ask it not to)
#08/21/2016 JMo updated dbCAN_hmm_file path and reworked mysql excute to load into feature_evidence based on feature_id which is contained in the sequence headers of the peptide files


use Getopt::Std;
use DBI;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

our $dbCAN_hmm_file = "/scripts/db/hmms/dbCAN/dbCAN-fam-HMMs.v5.txt";
our $DEBUG = 1;

my %arg;
&getopts('D:u:p:i:n:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $set_name = $arg{'n'};
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

sub run_dbCAN_hmmscan {
    my $pepfile = shift;
    system("hmmscan --cpu 10 $dbCAN_hmm_file $pepfile > $pepfile.dbCAN.hmmscan");
    return "$pepfile.dbCAN.hmmscan";
}

sub run_dbCAN_hmmsearch {
    my $pepfile = shift;
    system("hmmsearch --cpu 10 --domtblout $pepfile.dbCAN.domtblout --noali --acc -o $pepfile.dbCAN.out $dbCAN_hmm_file $pepfile");
    return "$pepfile.dbCAN.domtblout";
}

sub parse_hmmscan_output {
    my $hmmscan = shift;
    my @a;
    while (<>) {
	if(/^\/\//) { # when we find the end of a record
	    my $x = join("", @a); # join @a into one long string?
	    my ($q)=($x=~/^Query:\s+(\S+)/m); # look for the query acc
	    while ($x=~/^>> (\S+.*?\n\n)/msg) { # look for lines describing hit
		my $a = $&;
		my @c=split(/\n/,$a); # split this into individual lines
		$c[0]=~s/>> //; # grab the target accession
		for (my $i=3;$i<=$#c;$i++) { # skip the header lines
		    my @d=split(/\s+/,$c[$i]); # split the data line into fields
		    # here we're printing the prot acc, ev acc, i-Evalue, hmmfrom, hmmto, alifrom, alito 
		    print $q."\t".$c[0]."\t$d[6]\t$d[7]\t$d[8]\t$d[10]\t$d[11]\n" if $d[6]<1;
		}
	    }
	    @a=();
	} else {
	    push(@a,$_); # @a is all the lines in the record
	}
    }
    # sort by prot acc, alifrom, alito and feed into...
    my %b;
    while (<>) {
	chomp;
	my @a=split;
	next if $a[-1]==$a[-2]; # this is if alifrom and to are the same?
	push(@{$b{$a[0]}},$_); # make data structure %b with key of prot acc and value being an array of all hit info
    }
    foreach (sort keys %b) { # foreach prot acc
	my @a=@{$b{$_}};        # @a is the array of hits
	for(my $i=0;$i<$#a;$i++) { # foreach hit row
	    my @b=split(/\t/,$a[$i]);  # turn the info into fields
	    my @c=split(/\t/,$a[$i+1]); # turn the next info into fields, too
	    my $len1=$b[-1]-$b[-2]; # i think here we're seeing if the hits overlap
	    my $len2=$c[-1]-$c[-2];
	    my $len3=$b[-1]-$c[-2]; # this is b[high] - c[low], which, if >0 means an overlap
	    if ($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)) { # if the overlap is > 50% of the alilength
		if($b[2]<$c[2]) { # if b[evalue] < c[evalue]
		    splice(@a,$i+1,1); # remove c
		} else {
		    splice(@a,$i,1); # remove b
		}
		$i=$i-1; # since we removed an element, back up
	    }
	}
	foreach(@a) {
	    print $_."\n";
	}
    }
    # uniq
    open (IN,"all.hmm.ps.len");
    # read in this length file?
    my %b;
    while (<IN>) {
	chomp;
	my @a=split;
	$b{$a[0]}=$a[1]; # make a structure holding the hmm length
    }
    while (<>) {
	chomp;
	my @a=split; # split the hits rows into fields
	my $r=($a[4]-$a[3])/$b{$a[1]}; # calculate %of hmm in alignment
	print $_."\t".$r."\n"; # append this to the hit row.
    }
    while(<>) { # report hits with < 1e-5 if alilength >80, else 1e-3
	my @a=split(/\t/,$_);
	if (($a[-2]-$a[-3])>80) { 
	    print $_ if $a[2]<1e-5;
	} else {
	    print $_ if $a[2]<1e-3;
	}
    }
}

sub parse_hmmsearch_output {
    my $file = shift;
    open my $in, $file or die "Can't open $file: $!\n";
    my %dbcanhits;
    while (my $line = <$in>) {
	next if ($line =~ /^#/);
	chomp $line;
	my ($prot_acc,
	    undef,
	    $prot_len,
	    $ev_acc,
	    undef,
	    $ev_len,
	    $evalue,
	    $score,
	    $fullbias,
	    $domidx,
	    $totdoms,
	    $cevalue,
	    $ievalue,
	    $domscore,
	    $dombias,
	    $hmmfrom,
	    $hmmto,
	    $alifrom,
	    $alito,
	    $envfrom,
	    $envto,
	    $description) = split/\s+/, $line, 22;
	    # The dbCAN paper lists 3 criteria for a 'good' hit.
	    # length > 80 aa and e-value <1e-5
	    # length < 80 aa and e-value 
	    # HMM must be covered >50%
	push @{$dbcanhits{$prot_acc}}, [$ev_acc,
					$ievalue,
					$ev_len,
					$hmmfrom,
					$hmmto,
					$prot_len,
					$alifrom,
					$alito,
					($hmmto-$hmmfrom/$ev_len)] if (
	    ($ievalue < 1e-5 || ($ievalue < 1e-3 && ($alifrom-$alito+1) < 80)) &&
	    ($hmmto-$hmmfrom+1)/$ev_len >= 0.5);
    }
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
