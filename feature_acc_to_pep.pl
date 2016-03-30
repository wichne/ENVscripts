#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

&getopts('D:i:o:');
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

my @feat_accs;
if ($opt_i) {
    open my $in, $opt_i;
    while (my $line=<$in>) {
	my @f = split/\s+/, $line;
	push @feat_accs, $f[0];
    }
} else {
    while ($l = <STDIN>) {
	chomp $l;
	if (!$l) { last }
	$l =~ s/^>//;
	my @f = split/\s/, $l;
	push @feat_accs, $f[0];
    }
}

$fid_ref = &get_feature_id_by_accession($dbh, @feat_accs);
my @feat_ids;
foreach $ref (values %$fid_ref) {
    push @feat_ids, $ref->{'feature_id'};
}
my $protref = &get_protein_by_feature_id($dbh, @feat_ids);

foreach my $fid (keys %{$protref}) {
    my $acc = join "|", @{$protref->{$fid}->{'accessions'}};
    if (length($protref->{$fid}->{'product'}) == 0 ||
	!defined($protref->{$fid}->{'product'})) {
	warn "No product string for feature $fid ($acc). Skipping...";
	next;
    }
    my $desc;
    while (my ($k,$v) = each %{$protref->{$fid}->{'annotation'}}) {
	if ($k eq "product") { $desc = $v->[0] . " " . $desc; }
	else { 
	    foreach my $val (@$v) {
		$desc .= " [$k=$val]";
	    }
	}
    }
    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $protref->{$fid}->{'product'});
    $outfh->write_seq($seqo);
}
