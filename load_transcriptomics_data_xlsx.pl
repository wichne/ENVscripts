#!/usr/bin/perl
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;
use Spreadsheet::ParseXLSX;
use Spreadsheet::ParseExcel;
use strict;

$| = 1;

my %args;
&getopts('D:u:p:f:h', \%args);
if ($args{'h'}) {
    print "USAGE load_transcriptomics_data_xlsx.pl -f xlsxfile -D database -p dbpassword [ -u dbuser -h ]

-f path to the xlsx file that contains the data. Each sample should be on a separate tab with names that will sort in a useful manner (ie 0-fill your numbers). Also, there must be an accession or identifier by which the program can look up feature_id in the database in the '#Gene-ID' column.

-D the database to update
-u database user (optional) Default is current username.
-p the mysql password

-h This help message
";
    exit();
}
if (! $args{f}) { die "You need to specify a tsv transcriptomics infile with -f.\n"; }
my $expt = $args{f};
$expt =~ s/.*\/// if ($expt =~ /\//);
$expt =~ s/\..*//;

my $dbh = ENV::connect(\%args);

my $parser;
if ($args{f} =~ /.xlsx$/) {
    $parser = Spreadsheet::ParseXLSX->new;
} elsif ($args{f} =~ /.xls$/) {
    $parser = Spreadsheet::ParseExcel->new;
}
print STDERR "Parsing $args{f}...";
my $wkbk = $parser->parse($args{f});
print "Done\n";

my %DATA;
my %TOT;
my %LENGTH;

# create a datastructure that links feature_id to set_id
print STDERR "Getting genes...";
my %SET;
my $q = "SELECT s.set_id, m.feature_id FROM sequence_sets s, seq_set_link l, seq_feat_mappings m"
    . " WHERE l.seq_id=m.seq_id"
    . " AND s.set_id=l.set_id"
    . " AND s.is_current = 1";
my $r = $dbh->selectall_arrayref($q);
foreach my $row (@$r) {
    my ($set_id, $fid) = @$row;
    $SET{$fid} = $set_id;
}
print STDERR "done.\n";

# Let's look at the file, sheet by sheet
for my $wksht ($wkbk->worksheets() ) {
    my $sample = $wksht->get_name;
    print STDERR "Found worksheet $sample\n";
    if ($sample =~ /^Sheet\d+$/) {
	print STDERR "\tSkipping because it isn't named...\n";
	next;
    }
    
    my ( $row_min, $row_max ) = $wksht->row_range();
    my ( $col_min, $col_max ) = $wksht->col_range();
    
    # the first row is headers, or should be.
    # We're looking for:
    # #Gene-ID - holds gene accession
    # Chr      - holds bin
    # Length   - holds gene length in bp
    # Count    - holds read count
    my ($acc_col, $bin_col, $gene_length_col, $count_col);
    for my $col ($col_min .. $col_max) { 
	my $cell = $wksht->get_cell($row_min, $col);
	if (! defined $cell) { print STDERR "No cell: $row_min, $col\n"; next; }
	if ($cell->value =~ /\#Gene\-ID/i) {
	    $acc_col = $col;
	} elsif ($cell->value =~ /Chr/i) {
	    $bin_col = $col;
	} elsif ($cell->value =~ /Length/i) {
	    $gene_length_col = $col;
	} elsif ($cell->value =~ /Count/i) {
	    $count_col = $col;
	}
    }
    
    # Now parse through the data rows, pulling the cols of interest
    for my $row ( $row_min+1 .. $row_max ) {
	
	# grab the accession
	my $acc;
	my $fid;
	my $set_id;
	if (defined $wksht->get_cell($row, $acc_col) &&
	    $wksht->get_cell($row, $acc_col)->value) {
	    $acc = $wksht->get_cell($row, $acc_col)->value;
	}
	if ($acc =~ /^$args{D}\_(\d+)$/) {
	    $fid = $1;
	} elsif ($acc =~ /\|/) {
	    my @a = split/\|/, $acc;
	    foreach my $ac(@a) {
		if ($ac =~ /gnl\_/) {
		    $ac =~ s/gnl\_//;
		}
		$fid = get_feature_id_by_accession($dbh, $ac);
		if ($fid) { last }
	    }
	} else {
	    if ($acc =~ /gnl\_/) {
		$acc =~ s/gnl\_//;
	    }
	    $fid = get_feature_id_by_accession($dbh, $acc);
	}
	if (! $fid) {
	    warn "Couldn't find a feature_id for accession $acc. Skipping...\n";
	    next;
	}
	if (defined $SET{$fid}) {
	    $set_id = $SET{$fid};
	} else {
	    warn "Can't find a set_id for feature_id $fid. Skipping...\n";
	    next;
	}
	
	# skip non-CDS rows
	if (!(is_CDS($dbh, $fid))) {
	    warn "$acc is not a CDS, so we'll skip it...\n";
	    next;
	}
	
	# get the length
	if (defined $wksht->get_cell($row, $gene_length_col) && $wksht->get_cell($row, $gene_length_col)) {
	    $LENGTH{$fid} = $wksht->get_cell($row, $gene_length_col)->value;
	} else {
	    warn "No length for $acc?\n";
	    next;
	}
	
	# Get the count
	if (defined $wksht->get_cell($row, $count_col) && $wksht->get_cell($row, $count_col)) {
	    $DATA{$fid}->{$expt}->{$sample} = $wksht->get_cell($row, $count_col)->value;
	    $TOT{$set_id}->{$expt}->{$sample} += $DATA{$fid}->{$expt}->{$sample};
	} else {
	    $DATA{$fid}->{$expt}->{$sample} = 0;
	}
    }
}

open (my $set_count, ">${expt}_set_count.txt");
foreach my $fid (sort {$SET{$a} <=> $SET{$b} || $a<=>$b } keys %DATA) {
    print $set_count "$SET{$fid}\t";
    foreach my $ex (keys %{$DATA{$fid}}) {
	print $set_count "$ex\t";
	foreach my $samp (sort keys %{$DATA{$fid}->{$ex}}) {
	    print $set_count "$samp\t";
	    my $rel_count = 0;
	    if ($TOT{$SET{$fid}}->{$ex}->{$samp} > 0) {
		$rel_count = sprintf "%.2f", ($DATA{$fid}->{$ex}->{$samp}*100/$TOT{$SET{$fid}}->{$ex}->{$samp});
	    }
	    my $RPKM;
	    if (! (defined $DATA{$fid}->{$ex}->{$samp} ||
		   defined $LENGTH{$fid} ||
		   defined $TOT{$SET{$fid}}->{$ex}->{$samp}) ||
		$LENGTH{$fid} == 0 ||
		$TOT{$SET{$fid}}->{$ex}->{$samp} == 0) {
		$RPKM = 0;
	    } else {
		$RPKM = $DATA{$fid}->{$ex}->{$samp}/($LENGTH{$fid}/1000)/($TOT{$SET{$fid}}->{$ex}->{$samp}/1000000);
	    }
	    print $set_count $TOT{$SET{$fid}}->{$ex}->{$samp} . "\n";
#	    print "$SET{$fid}\t$fid\t$ex\t$samp\t" 
#		. $DATA{$fid}->{$ex}->{$samp} 
#	    . "\t" . $TOT{$SET{$fid}}->{$ex}->{$samp}
#	    . "\t$rel_count\t$RPKM\n";
	    &insert_transcriptomics_data($dbh, $fid, $ex, $samp, $DATA{$fid}->{$ex}->{$samp}, $RPKM);
##	    &add_abs_value($dbh, $fid, $expt, $samp, $DATA{$fid}->{$ex}->{$samp});
	}
    }
}

sub insert_transcriptomics_data {
    my $dbh = shift;
    my ($fid, $expt, $samp, $abs_val, $norm_val) = @_;
    
    my $q = "INSERT INTO transcriptomics_data (feature_id, experiment, sample, abs_value, norm_value) VALUES (?, ?, ?, ?, ?)";
    
    my $ret = $dbh->do($q, {}, ($fid, $expt, $samp, $abs_val, $norm_val));
    if (! defined $ret) { print STDERR "oops on $fid, $expt, $samp, $abs_val, $norm_val: $q\n"; }
}

sub add_abs_value {
    my $dbh = shift;
    my ($fid, $expt, $samp, $abs_val) = @_;

    my $q = "UPDATE transcriptomics_data set abs_value=? where feature_id=? and experiment=? and sample=?";
    
    my $ret = $dbh->do($q, {}, ($abs_val, $fid, $expt, $samp));
    if (! defined $ret) { print STDERR "oops2 on $fid, $expt, $samp, $abs_val : $q\n"; }
}


sub is_CDS {
    my $dbh = shift;
    my $fid = shift;
    my $q = "select feat_type from sequence_features where feature_id=$fid";
    my ($ft) = $dbh->selectrow_array($q);
    return $ft eq "CDS" ? 1 : 0;
}
