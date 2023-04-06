#!/usr/bin/perl

use DBI;
use Getopt::Std;
$| = 1;
my %args;
&getopts('D:i:n:o:',\%args);
my $db = $args{'D'};
my $setid = $args{'i'};
my $setname = $args{'n'};
my $outfile = $args{'o'};
my $host=$ENV{DBSERVER};

my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'access');

if ($setname) {
    $set_id_q = "SELECT set_id, name FROM sequence_sets WHERE name='$setname'";
    ($setid) = $dbh->selectrow_array($set_id_q);
} elsif ($setid) {
    $set_name_q = "SELECT set_name from sequence_sets where set_id=$setid";
    ($setname) = $dbh->selectrow_array($seq_name_q);
}

if ($outfile) {
    open OUT, ">$outfile" or die "Can't open $outfile for writing: $!\n";
    select OUT;
}

my @facc_src = ("IMG", "RAST", "PNNL");

my %DATA;
my $faccq = "SELECT f.feature_id, m.seq_id, feat_min, feat_max, strand, product"
    . " FROM sequence_features f, seq_feat_mappings m, seq_set_link l"
    . " WHERE l.set_id=$setid"
    . " AND m.seq_id=l.seq_id"
    . " AND f.feature_id=m.feature_id";
my $sth = $dbh->prepare($faccq);
$sth->execute or die "$faccq\n";
while (my $ref = $sth->fetchrow_hashref) {
    $DATA{$ref->{'feature_id'}}->{'location'} = $ref;
}

my $faccq = "SELECT x.feature_id, source, prefix, accession FROM feature_accessions x, seq_feat_mappings m, seq_set_link l"
    . " WHERE l.set_id=$setid"
    . " AND m.seq_id=l.seq_id"
    . " AND x.feature_id=m.feature_id";

$faccr = $dbh->selectall_arrayref($faccq);
foreach my $row (@$faccr) {
    push @{$DATA{$row->[0]}->{'accession'}->{$row->[1]}}, $row->[3];
}

# e.* is feature_id, feat_min, feat_max, program, ev_type, ev_accession, ev_min, ev_max, ev_length, score
my $feature_q = "SELECT e.* FROM feature_evidence e, seq_feat_mappings m, seq_set_link l"
    . " WHERE ev_type in ('TMH', 'SP', 'OMP')"
    . " AND l.set_id=$setid"
    . " AND m.seq_id=l.seq_id"
    . " AND e.feature_id=m.feature_id"
    . " ORDER BY feat_min";
my $sth = $dbh->prepare($feature_q);
$sth->execute;
while (my $ref = $sth->fetchrow_hashref) {
    push @{$DATA{$ref->{'feature_id'}}->{'evidence'}->{$ref->{'ev_type'}}->{$ref->{'ev_accession'}}}, $ref;
}

my $ann_q = "SELECT a.feature_id, i.qualifier, a.value, a.source, a.rank"
    . " FROM feature_annotations a, INSDC.qualifier i, seq_feat_mappings m, seq_set_link l"
    . " WHERE l.set_id=$setid"
    . " AND m.seq_id=l.seq_id"
    . " AND a.feature_id=m.feature_id"
    . " AND i.id=a.data_type_id"
    . " ORDER BY rank DESC";
my $sth = $dbh->prepare($ann_q);
$sth->execute;
while (my $ref = $sth->fetchrow_hashref) { 
    push @{$DATA{$ref->{'feature_id'}}->{'annotation'}->{$ref->{'qualifier'}}->[$rank]}, $ref;
}

print "Feature_id\tIMG acc\tRAST acc\tPNNL acc\tSeq_id\tMin\tMax\tStrand\tProduct\tGene_sym\tEC#\tTMH count\tTMH loc\tLepB SP score\tLepB SP\tLspA SP score\tLspA SP\tTAT SP score\tTAT SP\tOMP sig score\tOMP sig\n";
foreach my $feature_id (sort keys %DATA) {
    my $ref = $DATA{$feature_id};
    my @data = ($feature_id);
    foreach my $at (@facc_src) {
	my @T;
	foreach my $r (@{$ref->{'accession'}->{$at}}) {
	    push @T, $r;
	}
	push @data, join(" ", @T);
    }
    my $strand = $ref->{location}->{strand} == -1 ? "-" : $ref->{location}->{strand} == 1 ? "+" : ".";
    push @data, ($ref->{location}->{seq_id},
		 $ref->{location}->{feat_min},
		 $ref->{location}->{feat_max},
		 $strand
    );
    
    my $prod;
    foreach my $r (@{$ref->{'annotation'}->{'product'}}) {
	foreach $t(@$r) {
	    if ($t->{source}) {
		$prod .= "$t->{source}=$t->{value};"
	    } else { $prod .= $t->{value} }
	}
    }
    my $gene;
    foreach my $r (@{$ref->{'annotation'}->{'gene'}}) {
	foreach $t(@$r) {
	    $prod .= "$t->{source}=$t->{value};"
	}
    }
    my $ec;
    foreach my $r (@{$ref->{'annotation'}->{'EC_number'}}) {
	foreach $t(@$r) {
	    $prod .= "$t->{source}=$t->{value};"
	}
    }
    push @data, ($prod, $gene, $ec);

    my @tm;
    push @data, scalar(@{$ref->{'evidence'}->{'TMH'}->{'TMhelix'}});
    foreach my $t (@{$ref->{'evidence'}->{'TMH'}->{'TMhelix'}}) {
	push @tm, $t->{'feat_min'} . ".." . $t->{'feat_max'};
    }
    push @data, join(" ", @tm);

    push @data, ($ref->{'evidence'}->{'SP'}->{'SpI'}->[0]->{'score'},
		 substr($ref->{'location'}->{'product'}, 0, $ref->{'evidence'}->{'SP'}->{'SpI'}->[0]->{'feat_max'}));
    push @data, ($ref->{'evidence'}->{'SP'}->{'SpII'}->[0]->{'score'},
		 substr($ref->{'location'}->{'product'}, 0, $ref->{'evidence'}->{'SP'}->{'SpII'}->[0]->{'feat_max'}));
    push @data, ($ref->{'evidence'}->{'SP'}->{'Tat_signal_seq'}->[0]->{'score'},
		 substr($ref->{'location'}->{'product'}, 0, $ref->{'evidence'}->{'SP'}->{'Tat_signal_seq'}->[0]->{'feat_max'}));

    my @omp;
    foreach my $oacc (keys %{$ref->{'evidence'}->{'OMP'}}) {
	push @omp, ($oacc . "::" . $ref->{'evidence'}->{'OMP'}->{$oacc}->[0]->{'feat_min'} . ".." . $ref->{'evidence'}->{'OMP'}->{$oacc}->[0]->{'feat_max'}
#		    substr($ref->{'location'}->{'product'},
#			   $ref->{'evidence'}->{'OMP'}->{$oacc}->[0]->{'feat_min'},
#			   $ref->{'evidence'}->{'OMP'}->{$oacc}->[0]->{'feat_max'})
	);
    }
    push @data, join(" ", @omp);

    print join("\t", @data) . "\n";
}

