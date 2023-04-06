#!/usr/bin/perl
use lib "/home/wcnelson/devel/";
use ENV;
use DBI;
use Bio::SeqFeature::Generic;
use Carp;
use strict;
use Getopt::Std;

my %args;
&getopts('D:I:i:u:p:j:h',\%args);

if ($args{h}) {
    print "";
}

my @ids;
if ($args{i}) {
    @ids = $args{i};
} elsif ($args{I}) {
    open my $ids, $args{I};
    while (<$ids>) {
	push @ids, $_;
    }
} else {
    @ids = <STDIN>;
}

my $dbh = DBI->connect("dbi:mysql:server=localhost;database=$args{D}", $args{u}, $args{p});
my $joint_q = "SELECT COUNT(sf.feature_id)"
    . " FROM seq_feat_mappings sfm, sequence_features sf, sequence_accessions sa"
    . " WHERE sa.seq_accession=?"
    . " AND sfm.seq_id=sa.seq_id"
    . " AND sfm.feature_id=sf.feature_id"
    . " AND sf.feat_type_id=12"; # feat_type_id 12 is 'scaffold_joint'
my $joint_sth = $dbh->prepare($joint_q);

my $seq_q = "SELECT sequence"
    . " FROM sequences s, sequence_accessions sa"
    . " WHERE sa.seq_accession=?"
    . " AND sa.seq_id = s.seq_id";
my $seq_sth = $dbh->prepare($seq_q);

my $regex = "N{20,}";
if ($args{j}) {
    $regex = $args{j};
}

# examine scaffolds one at a time
foreach my $id (@ids) {
    chomp $id;
    # check to see if scaffold already has scaffold_joint features
    my $jrv = $joint_sth->execute($id);
    if (! defined $jrv) { croak $dbh->errstr }
    my ($joint_count) = $joint_sth->fetchrow_array;

    if ($joint_count > 0) {
	carp "Scaffold with seq_accession $id already has scaffold_joint rows";
	next;
    }

# get the scaffold sequence
    my $srv = $seq_sth->execute($id);
    my ($seq) = $seq_sth->fetchrow_array;

# look for scaffold joint sequence
# default is string of 20 or more N's
    while ($seq =~ /$regex/gi) {
	my $lo_coord = length($`) + 1;
	my $hi_coord = pos $seq;
# load these coords as scaffold_joint features
	my $feato = new Bio::SeqFeature::Generic(-start   => $lo_coord,
						 -end     => $hi_coord,
						 -strand  => 0,
						 );
	print "$id\t$lo_coord\t$hi_coord\n";
	&load_SeqFeature($dbh, $id, 12, $feato, 0000730);
    }
}
