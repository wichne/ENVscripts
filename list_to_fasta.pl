#!/usr/bin/perl

#04/12/2016 WCN Based on set_to_seq.pl. Provide list of feature_ids. Produce fasta file.

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper qw(Dumper);

&getopts('D:U:P:i:I:pso:h');
if ($opt_h) {
    die "USAGE:
list_to_fasta.pl -D <db> -P <dbpswd> [ -U <dbusername> ] [ -i <feature_id> OR -I <file_of_feature_ids> ] [ -p OR -s ] [-o <output_file> ] [ -h ]\n";
}

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $user = $opt_U ? $opt_U : $ENV{USER};
my $pswd = $opt_P || die "Must provide password with -P";
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$opt_D", $user, $pswd);

my @feature_ids;
if ($opt_i) {
    @feature_ids = $opt_i;
} elsif ($opt_I) {
    open my $ids, $opt_I;
    my $data_errors;
    while (<$ids>) {
	chomp;
	if ($_ =~ /[\D]/) { $data_errors .= "$_; " }
	push @feature_ids, $_;
    }
    if ($data_errors) { die "Invalid data in $opt_I: $data_errors\n"; }
}

my $outfh;
if ($opt_o) {
    $outfh = Bio::SeqIO->new(-file => ">$opt_o",
			     -format => 'fasta');
} else {
    $outfh = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => 'fasta');
}

my $type = $opt_s ? "sequence" : "protein";

if ($type eq "sequence") {
    my $seq_ref = get_sequence_by_feature_id($dbh, @feature_ids);
    foreach my $feature_id(keys %$seq_ref) {
	my $acc = get_locus_tag_by_feature_id($dbh, $feature_id)->[0];
	my $set = seq_id_to_set_name($dbh, $seq_ref->{$feature_id}->{'seq_id'});
	my $min = $seq_ref->{$feature_id}->{'min_partial'} ? "<" . $seq_ref->{$feature_id}->{'feat_min'} : $seq_ref->{$feature_id}->{'feat_min'};
	my $max = $seq_ref->{$feature_id}->{'max_partial'} ? ">" . $seq_ref->{$feature_id}->{'feat_max'} : $seq_ref->{$feature_id}->{'feat_max'};
	my $coords = $min . ".." . $max;
	if ($seq_ref->{$feature_id}->{'strand'} == -1) { $coords = "complement(" . $coords . ")"}
	my $desc = $set  . " " . $coords . " " . $seq_ref->{$feature_id}->{'annotation'}->{'product'};
	my $seq = $seq_ref->{$feature_id}->{'gene_seq'};
	my $seqo = Bio::Seq->new(-display_id => $acc,
				 -desc => $desc,
				 -seq => $seq);
	$outfh->write_seq($seqo);
    }
} elsif ($type eq "protein") {


foreach my $fid (@feature_ids) {
    get_sequence_by_feature_id
	get_features_by_feature_id
	
#    my $acc = join("|", values %{$protref->{$fid}->{'accessions'}});
    my $locus = $protref->{$fid}->{'accessions'}->{'locus_tag'};
    my $db = $opt_D;
    my $acc = "$db|$fid|$locus";
    
    my ($seq_id, $min, $max, $strand) = ($protref->{$fid}->{'location'}->{'seq_id'},
					 $protref->{$fid}->{'location'}->{'feat_min'},
					 $protref->{$fid}->{'location'}->{'feat_max'},
					 $protref->{$fid}->{'location'}->{'strand'});
    if (!($min && $max) || ($min >= $max)) { warn "Why $min/$max for $fid, $acc?\n";}
    my $seq = &get_subseq($seq_id, $min, $max, $strand);
    
#     if (length($protref->{$fid}->{'product'}) == 0 ||
# 	!defined($protref->{$fid}->{'product'})) {
# 	warn "No product string for feature $fid ($acc). Skipping...";
# 	next;
#     }
    my $setdesc = &get_set_desc($dbh, $setid);
    my $desc = "[$setdesc]";
#    while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
#	my ($k) = sort {$a<=>$b} keys %$vref;
#	$desc .= "$qual=$vref->{$k}->[0]->{value};";
#}	

    my $seqo = Bio::Seq->new(-display_id => $acc,
			     -desc => $desc,
			     -seq => $seq);
    $outfh->write_seq($seqo);
}
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

sub get_subseq {
    my ($seq_id, $min, $max, $strand) = @_;
    my $subseq = $SEQ->{$seq_id}->trunc($min, $max);
    if ($strand < 0 || $strand eq "-") {
	return $subseq->revcom->seq;
    } else {
	return $subseq->seq;
    }
}
