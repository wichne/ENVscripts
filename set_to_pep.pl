#!/usr/bin/perl
#Edited 05/17/2015 by Jen in order to incorporate only the databasename and feature id and organism description
#in sequence headers

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

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

foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
#    my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
    my $db = $opt_D;
    my $acc = "$db|$fid"; 

    if (length($protref->{$fid}->{'product'}) == 0 ||
	!defined($protref->{$fid}->{'product'})) {
	warn "No product string for feature $fid ($acc). Skipping...";
	next;
    }
    my $setdesc = &get_set_desc($dbh,$setid);
    my $desc = "[$setdesc]";
#    while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
#	my ($k) = sort {$a<=>$b} keys %$vref;
#	$desc .= "$qual=$vref->{$k}->[0]->{value};";
#}	

    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $protref->{$fid}->{'product'});
    $outfh->write_seq($seqo);
}

sub get_set_desc {
    my $dbh = shift;
    my $setid = shift;
    my $ssq = "SELECT description FROM sequence_sets WHERE set_id=$setid";
    return $dbh->selectcol_arrayref($ssq)->[0];
}

sub get_mainroles {
    my $dbh = shift;
    my $mrq = "select * from egad.mainrole";
    return $dbh->selectall_hashref($mrq, 'id');
}

sub get_subroles {
    my $dbh = shift;
    my $mrq = "select * from egad.subrole";
    return $dbh->selectall_hashref($mrq, 'id');
}
