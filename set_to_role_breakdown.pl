#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use strict;

our($opt_D, $opt_i, $opt_n, $opt_o);

&getopts('D:i:n:o:');
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'access');

my $setid = $opt_i;
my $setname = $opt_n;
my $outfh;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");

my $MR = &get_mainroles($dbh);
my $SR = &get_subroles($dbh);

my %MAINROLE;
my %SUBROLE;
foreach my $fid (keys %$protref) {
    my ($mainrole, $subrole);
    my $q = "SELECT subrole, LEFT(subrole,2) AS mainrole FROM egad.ev_role_link erl, feature_evidence fe"
	. " WHERE fe.feature_id=$fid"
	. " AND erl.ev_acc=fe.ev_accession";
    my $rret = $dbh->selectall_hashref($q, 'subrole');
    if (scalar(keys %$rret)> 1) {
	print STDERR "$fid $protref->{$fid}->{'accessions'}->[0] evidence points to multiple roles: " . join(";",(keys %$rret)) . "\n";
    } elsif (! %$rret) {
	$MAINROLE{99}++;
	$SUBROLE{9999}++;
    }
    $MAINROLE{$rret->{[keys %$rret]->[0]}->{'mainrole'}}++;
    $SUBROLE{$rret->{[keys %$rret]->[0]}->{'subrole'}}++;
    if (defined $protref->{$fid}->{'accessions'}) {
	my $acc = join "|", @{$protref->{$fid}->{'accessions'}};
#	my $acc = {$protref->{$fid}->{'accessions'}};
	print join("\t",($acc,[keys %$rret]->[0],$MR->{$rret->{[keys %$rret]->[0]}->{'mainrole'}}->{'name'},$SR->{$rret->{[keys %$rret]->[0]}->{'subrole'}}->{'name'})) . "\n";
    }
}

# Print the output
# First the mainrole breakdown
# Print a header line
print "Bin\tCDS\t";
my @mrs;
foreach my $mrid (sort {$a cmp $b} keys %$MR) {
    push @mrs, $MR->{$mrid}->{'name'};
}
print join("\t", @mrs), "\n";

# Now print the data
my $cdscount = scalar(keys%$protref);
print "$setname\t$cdscount\t";
my @counts;
foreach my $mrid (sort {$a cmp $b} keys %$MR) {
    push @counts, $MAINROLE{$mrid}/$cdscount;
}
print join("\t", @counts),"\n";

## Now the subrole breakdown
# Print a header line
my (@smrs, @srs);
foreach my $srid (sort {$a cmp $b} keys %$SR) {
    $srid =~ /(.{2})/;
    my $mrid = $1;
    push @smrs, $MR->{$mrid}->{'name'};
    push @srs, $SR->{$srid}->{'name'};
}
print "\t\t";
print join("\t", @smrs), "\n";
print "Bin\tCDS\t";
print join("\t", @srs), "\n";
# Now print the data
print "$setname\t$cdscount\t";
my @counts;
foreach my $srid (sort {$a cmp $b} keys %$SR) {
    push @counts, $SUBROLE{$srid}/$cdscount;
}
print join("\t", @counts),"\n";
