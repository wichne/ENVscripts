#!/usr/bin/env perl

use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

&getopts('D:i:n:o:O');
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", 'access', 'mySQL@cce55');

my $setid = $opt_i;
my $setname = $opt_n;
#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }
if ($setid && !$setname) { $setname = &set_id_to_set_name($dbh, $setid); }

my $seqref = &get_sequences_by_set_id($dbh, $setid);

my $outfh;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
} elsif ($opt_O) {
    $outfh = Bio::SeqIO->new(-file => ">${setname}.fasta",
			     -format => 'fasta');
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
}

foreach my $seqid (keys %$seqref) {
    my $acc = join "|", sort @{$seqref->{$seqid}->{'accessions'}};
    
    while (my ($k,$v) = each %{$protref->{$fid}->{'annotation'}}) {
	foreach my $val (@$v) {
	    $desc .= "$k=$val;";
	}
    }

    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $seqref->{$seqid}->{'sequence'});
    $outfh->write_seq($seqo);
}

exit();
