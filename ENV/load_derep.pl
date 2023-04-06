#!/usr/bin/perl -w

###
# Input
#    fasta files (or headers) from prodigal, with seq_accession that exists in db already, start, end, direction/frame
# Example prodigal header:
# >ITZX_contig-100_0_1 # 2 # 772 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.554
#  -----------------    --- ----  ---         ----------
#     contig_acc       start end  dir          partial
#
#
#
#    usearch -derep-fulllength .uc file showing clusters of identical proteins
# Example:
# S       450     254     *       *       *       *       *       HLSNC17_01410   *
# H       450     254     100.0   *       0       0       *       ITZX_contig-100_302_9   HLSNC17_01410
# H       450     254     100.0   *       0       0       *       ITZ_contig-100_609_22   HLSNC17_01410
# H       450     254     100.0   *       0       0       *       ITZY_contig-100_71_14   HLSNC17_01410
# H       450     254     100.0   *       0       0       *       fig|6666666.44075.peg.281       HLSNC17_01410
# H       450     254     100.0   *       0       0       *       fig|6666666.44075.peg.3539      HLSNC17_01410

# The first line (S) lists the representative accession
# The following H lines list all the members
# For each set, a feature will be loaded
# The accessions listed in field [8] will be loaded in the feature_accessions table
# If the accession is associated with a contig/pos from the first input file, seq_feat_mappings is also loaded
# Should also add in a check to see if any of the accessions in [8] are in the db already, and pull the feature_id
# from that existing record; also check to see if multiple accessions yield multiple feature_ids - that would be an issue.

# If the set does not include an accession linked to a contig/pos, dump to a file of features to be added by another mechanism

# Keep track of numbers and types of insertions happening.

use DBI;
use Getopt::Long;
use Bio::SeqIO;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

my ($db, $user, $password, @fastafiles, $derep_file);
&GetOptions('db=s' => \$db,
	    'user=s' => \$user,
	    'password=s' => \$password,
	    'fasta=s' => \@fastafiles,
	    'input=s' => \$derep_file);

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $logfile = "load_derep.$$.log";
open LOG, ">$logfile" or warn "Can't write to logfile $logfile:! \n";

@fastafiles = split(/,/,join(',',@fastafiles));

# Read in the coordinate mapping information
my %MAP;

foreach my $file (@fastafiles) {
    print LOG "$file...\n";
    my $sfo = Bio::SeqIO->new(-file=>$file);
    while (my $seqo = $sfo->next_seq) {
	my ($strand, $phase);
	my (undef, $low, $high, $frame, $desc) = split /\s?\#\s/, $seqo->desc;
	($phase, $strand) = $frame < 0 ? (abs($frame), -1) : ($frame, 1);
	$phase = $phase > 3 ? $phase - 4 : $phase - 1;
	$phase = 0 if ($phase < 0);
	$seqo->display_id =~ /((.*contig\-\d+\_\d+)\_\d+)/;
	my ($contig_acc, $feat_acc) = ($2, $1);
	$desc =~ /partial\=(\d)(\d)/;
	my ($min_partial, $max_partial) = ($1, $2);
	$MAP{$feat_acc} = {'seq_accession' => $contig_acc,
			   'feat_min' => $low,
			   'feat_max' => $high,
			   'strand' => $frame,
			   'phase' => $phase,
			   'min_partial' => $min_partial,
			   'max_partial' => $max_partial,
			   'product' => $seqo->seq};
    }
    print scalar(keys %MAP) . "\n";
}

# Read in the derep clusters
my %FEAT;
open my $in, $derep_file or die "Can't open $derep_file: $!\n";
while (my $line = <$in>) {
    chomp $line;
    my @f = split/\t/, $line;
    if ($f[0] eq "S") {
	push @{$FEAT{$f[8]}}, $f[8];
    } elsif ($f[0] eq "H") {
	if(! defined $FEAT{$f[9]}) { warn "No feat instantiated for $f[9]!\n";}
	push @{$FEAT{$f[9]}}, $f[8];
    }
}
close $in;

# Now for the meat: load the database
my $feat_insert;
my $acc_insert;
my $map_insert;
my $existing_feats;
my $existing_accs;
my $nomap;

my $sequence_features_i = "INSERT sequence_features"
    . " (feat_type_id, feat_type, product, inserted_by, date_inserted)"
    . " VALUES (6, 'CDS', ?, USER(), NOW())";
my $sth1 = $dbh->prepare($sequence_features_i);
my $feat_acc_i = "INSERT feature_accessions"
    . " (feature_id, source, prefix, accession)"
    . " VALUES (?, ?, ?, ?)" ;
my $sth2 = $dbh->prepare($feat_acc_i);
my $seq_feat_mappings_i = "INSERT seq_feat_mappings"
    . " (seq_id, feature_id, feat_min, feat_max, strand, phase, min_partial, max_partial)"
    . " VALUES (?,?,?,?,?,?,?,?)";
my $sth3 = $dbh->prepare($seq_feat_mappings_i);
my $feat_acc_q = "SELECT feature_id FROM feature_accessions"
    . " WHERE accession=? AND prefix=? AND source=?";
my $sth4 = $dbh->prepare($feat_acc_q);

foreach my $r (keys %FEAT) {
    # We need to first see if any of the accessions are already present and
    # associated with a feature_id
    my %ACC;
    my %FID;
    my $product;
    foreach my $racc (@{$FEAT{$r}}) {
	if (defined $MAP{$racc}->{product} && ! $product) { 
	    $product = $MAP{$racc}->{product};
	}
	if ($racc =~ /\|/) {
	    $racc =~ /(.*)\|(.*)/;
	    ($ACC{$racc}->{prefix}, $ACC{$racc}->{accession}) = ($1, $2);
	} else {
	    $ACC{$racc}->{prefix} = "";
	    $ACC{$racc}->{accession} = $racc;
	}
	$ACC{$racc}->{source} = $ACC{$racc}->{prefix} eq "fig" ? "RAST"
	    : $ACC{$racc}->{accession} =~ /^HL.*\_\d+/ ? "IMG"
	    : "PNNL";
	$sth4->execute($ACC{$racc}->{accession}, $ACC{$racc}->{prefix}, $ACC{$racc}->{source});
	while (my ($feat_id) = $sth4->fetchrow_array) {
	    push @{$ACC{$racc}->{feature_id}}, $feat_id;
	    $FID{$feat_id}++;
	}
    }

    my $feat_id;
    if (scalar(keys %FID) == 1) {
	($feat_id) = keys %FID;
	$existing_feats++;
    } elsif (! %FID) {
	# insert the feature
	my $row1 = $sth1->execute($product);
	$feat_id = $dbh->last_insert_id("%", "%", "", "");
	if (! $feat_id) { die "Problem with seq_feature insert at $r\n";}
	$feat_insert++;
    } else {
	print LOG "Multiple feature ids assigned to unique feature.";
	print LOG join "\t",@{$FEAT{$r}};
	print LOG join "\t", keys %FID;
	warn "Multiple feature ids assigned to unique feature.";
	warn join "\t",@{$FEAT{$r}};
	warn join "\t", keys %FID;
	die;
    }

    my $mapped = 0;
    foreach my $kk (keys %ACC) {
	if (! defined $ACC{$kk}->{feature_id}) {
	    $sth2->execute($feat_id, $ACC{$kk}->{source}, $ACC{$kk}->{prefix}, $ACC{$kk}->{accession});
	    $acc_insert++;

	    if (defined $MAP{$kk}->{feat_min}) { # why is defined test failing?
		# map the feature
		my $seq_id = get_seq_id_by_accession($dbh, $MAP{$kk}->{seq_accession});
		$sth3->execute($seq_id, $feat_id,
			       $MAP{$kk}->{feat_min},
			       $MAP{$kk}->{feat_max},
			       $MAP{$kk}->{strand},
			       $MAP{$kk}->{phase},
			       $MAP{$kk}->{min_partial},
			       $MAP{$kk}->{max_partial});
		$mapped++;
		$map_insert++;
	    }
	} else { 

	    $existing_accs++;
	}
    }
    if (! $mapped) {
	print LOG "NOMAP: " . join(",", keys %ACC) . "\n";
	$nomap++;
    }
}

print LOG "Inserted $feat_insert features, $acc_insert accessions, and $map_insert mappings\n";
print LOG " $nomap features didn't map to existing contigs\n";
print LOG "There were $existing_feats features and $existing_accs accessions already in the database\n";

exit;
