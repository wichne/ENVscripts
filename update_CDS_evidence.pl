#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

# read in a list of accs or feat_ids (or a set_id?)
# make a pep fasta file of all the features
# Run:
# hmmer v TIGRFAMs
# hmmer v Pfam
# hmmer v dbCan
# signalp
# tmhmm
# lipop

my %args;
getopts('D:f:u:p:', \%args);
my $db = $args{'D'};
my $pwd = $args{'p'};

my $dbh = &connect(\%args);

my @ids;
my $fh;
if ($args{'f'}) { 
    open(my $fh, $args{'f'});
} else {
    $fh = *STDIN;
}

while (my $line = <$fh>) {
	chomp $line;
	$line =~ s/^\s*//;
	@line = split(/\s+/, $line);
	push @ids, $line[0];
}
close $fh;

if (!@ids) { die "No ids provided by -f or STDIN.\n"; }

# get the feature_ids
my @fids;
foreach my $id (@ids) {
    my $fid;
    if ($id =~ /^\d+$/) { $fid = $id }
    else { $fid = get_feature_id_by_accession($dbh, $id) }
    push @fids, $fid;
}

# build the fasta file
print STDERR "Building the fasta file...\n";
my $tmp_fasta = "$$.tmp.fasta";
my $outfh = Bio::SeqIO->new(-file => ">$tmp_fasta", -format => 'fasta');
my $features = get_features_by_feature_id($dbh, @fids);
foreach my $fid (keys %$features) {
    my $seqo = Bio::Seq->new(-display_id => $fid,
			     -seq => $features->{$fid}->{'product'});
    $outfh->write_seq($seqo);
}

$hmmer = "hmmsearch --acc --noali --cpu 16";
$TIGRFAM_db = "/scripts/db/hmms/TIGRFAM/TIGRFAMs_14.0_HMM.LIB";
$Pfam_db = "/scripts/db/hmms/Pfam/Pfam-A.hmm";
$dbCAN_db = "/scripts/db/hmms/dbCAN/dbCAN.hmm";

print STDERR "Running TIGRFAMs...\n";
system($hmmer . " -o $$.TIGRFAM.out --domtblout $$.TIGRFAM.domtblout --cut_nc $TIGRFAM_db $tmp_fasta");
system("load_hmms_to_ENV.pl -D $db -P $pwd -i $$.TIGRFAM.domtblout") or warn "$!\n";
# unlink "$$.TIGRFAM.out", "$$.TIGRFAM.domtblout";

print STDERR "Running Pfams...\n";
system($hmmer . " -o $$.Pfam.out --domtblout $$.Pfam.domtblout --cut_nc $Pfam_db $tmp_fasta");
system("load_hmms_to_ENV.pl -D $db -P $pwd -i $$.Pfam.domtblout") or warn "$!\n";
# unlink "$$.Pfam.out", "$$.Pfam.domtblout";

print STDERR "Running dbCAN...\n";
system($hmmer . " -o $$.dbCAN.out --domtblout $$.dbCAN.domtblout --cut_nc $dbCAN_db $tmp_fasta");
system("load_hmms_to_ENV.pl -D $db -P $pwd -i $$.dbCAN.domtblout") or warn "$!\n";
# unlink "$$.dbCAN.out", "$$.dbCAN.domtblout";

print STDERR "Running tmhmm...\n";
my $tmhmm = "tmhmm";
system("cat $tmp_fasta | $tmhmm > $$.tmhmm.out");
system("load_tmhmm.pl -D $db -P $pwd -i $$.tmhmm.out") or warn "$!\n";
# unlink "$$.tmmhmm.out";

print STDERR "Running lipop...\n";
my $lipop = "LipoP";
system("cat $tmp_fasta | $lipop > $$.lipop.out");
system("load_lipop.pl -D $db -P $pwd -i $$.lipop.out") or warn "$!\n";
# unlink "$$.lipop.out";

print STDERR "Running signalp...\n";
my $signalp = "signalp -t gram- -c 100";
system("cat $tmp_fasta | $signalp > $$.signalp.out");
system("load_signalp.pl -D $db -p $pwd -i $$.signalp.out") or warn "$!\n";
# unlink "$$.signalp.out";

#unlink $tmp_fasta);

