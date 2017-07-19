#!/usr/bin/perl
use lib $ENV{'ENVLIB'};
use ENV;
use Getopt::Std;
use strict;

my %args;
&getopts('D:p:i:d', \%args);
my $file = $args{'i'} || die "provide input with -i\n";
open (my $map, $file);
my $DEBUG = $args{'d'};

my $log_file = "load_mummer_remap.$$.log";
open (LOG, ">$log_file");

my $dbh = connect(\%args);

# 0 - fasta file
# 1 - seq_acc
# 2 - feat_type
# 3 - locus_tag
# 4 - low-coord
# 5 - hi-coord
# 6 - strand
# 7 - ref fasta file
# 8 - ref seq_acc
# 9 - ref feat_type
# 10 - ref locus_tag
# 11 - ref low coord
# 12 - ref high coord
# 13 - ref strand
# 14 - problems

my %seq_ids;
while (my $line = <$map>) {
    chomp $line;
    my @f = split(/\t/, $line);
    if ($f[0] eq "toDataSource") { next } # this is the header line

    my $ref_fid = $f[10] ? &get_feature_id_by_accession($dbh, $f[10]) : undef;
    my $seq_id = &get_seq_id_by_seq_accession($dbh, $f[1]);
    $seq_ids{$seq_id}++;

    # if the start appears to have been edited for the reference, adjust for this.
    my ($start, $end) = ($f[4], $f[5]);
    my $pseudo = 0;
    if ($f[14] =~ /mapped qry start=(\d+)/) {
	if ($f[6] >= 0 ) { $start = $1 } else { warn "for $f[1]::$f[10], check the end3\n"; }
    }
    if ($f[14] =~ /mapped qry end=(\d+)/) {
	if ($f[6] <=0 )  { $end = $1 } else {  warn "for $f[1]::$f[10], check the end3\n"; }
    }
    if (($end - $start + 1) % 3 != 0) { $pseudo = 5 }

    # Do we have a segmented gene?
    my $trans_string;
    if ($f[14] =~ /((complement\()?join\([^\)]+\)\)?)/) { $trans_string = $1 }
    else { $trans_string = $f[6] >= 0 ? "$start..$end" : "complement($start..$end)";}

    # The easy case - a perfect mapping
    if ($f[0] !~ /^#/) {
	if ($ref_fid) {
	    # insert a row linking the $f[10] feat_id to the mol/coords/strand in $f[1]/$f[4]-$f[5]/$f[6]
	    my $sfm_insert_i = "insert into seq_feat_mappings (seq_id, feature_id, feat_min, feat_max, strand, phase, min_partial, max_partial, pseudo, translation_coords)"
		. " VALUES ($seq_id, $ref_fid, $start, $end, '$f[6]', '0', 0, 0, '$pseudo', '$trans_string')";
	    print "$sfm_insert_i\n";
	    $dbh->do($sfm_insert_i) unless ($DEBUG);
	    print LOG "$f[10] ($ref_fid) mapped to $f[1] ($seq_id) $start-$end ($f[6])\n";
	} else {
	    # there is no old feature, so we must create a new one
	    my $seq_feat_i = "insert into sequence_features (feat_type, is_current)"
		." values ('$f[2]', 1)";
	    $dbh->do($seq_feat_i) unless ($DEBUG);
	    my $new_fid = $dbh->last_insert_id("%", "%", "", "");
	    my $sfm_insert_i = "insert into seq_feat_mappings (seq_id, feature_id, feat_min, feat_max, strand, phase, min_partial, max_partial, pseudo, translation_coords)"
		. " VALUES ($seq_id, $new_fid, $start, $end, '$f[6]', '0', 0, 0, '$pseudo', '$trans_string')";
	    print "$sfm_insert_i\n";
	    $dbh->do($sfm_insert_i) unless ($DEBUG);
	    print LOG "Inserted new feature ($f[2]) at $f[1] $start-$end\n";
	}
    } else {
        # rows with problem mappings start with "#"
        # if the row starts with a "#", then there is a significant difference in the 
	# gene model. The assumption is that the query genes will be more complete than
	# the ref genes, so insert a new feature.
	if (! $ref_fid) {
	    my $seq_feat_i = "insert into sequence_features (feat_type, is_current)"
	    ." values ('$f[2]', 1)";
	    $dbh->do($seq_feat_i) unless ($DEBUG);
	    $ref_fid = $dbh->last_insert_id("%", "%", "", "");
	}
	my $sfm_insert_i = "insert into seq_feat_mappings (seq_id, feature_id, feat_min, feat_max, strand, phase, min_partial, max_partial, pseudo, translation_coords)"
	    . " VALUES ($seq_id, $ref_fid, $start, $end, '$f[6]', '0', 0, 0, '$pseudo', '$trans_string')";
	    print "$sfm_insert_i\n";
	$dbh->do($sfm_insert_i) unless ($DEBUG);
	# if there's old accessions, link them to this new feature
#	if ($ref_fid) {
#	    my $acc_i = "insert into feature_accessions (feature_id, source, prefix, accession) select $new_fid, source, prefix, accession from feature_accessions where feature_id=$ref_fid";
#	    $dbh->do($acc_i);
#	}
#	print LOG "Inserted new feature $new_fid at $f[1] $start-$end to replace $f[10] ($ref_fid)\n";
	print LOG "Linked $f[10] ($ref_fid) to $seq_id at $start-$end\n";
    }
    # Need some code dealing with if the qry (new) features are loaded in the db
    # I guess they should be deprecated as the old features are mapped forward
}

# call rewrite seqs
my $rws_cmd = "rewrite_product.pl -db $args{D} -password $args{p}";
foreach my $sid (keys %seq_ids) {
#    print LOG "Running rewrite_product on $sid\n";
#    system($rws_cmd . " -seq_id $sid");
}

exit();

