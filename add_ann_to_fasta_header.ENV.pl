#!/usr/bin/perl

use DBI;
use Getopt::Std;
use strict;
$| = 1;
my %arg;
&getopts('D:u:p:f:', \%arg);

my $file = $arg{'f'};
open my $in, $file or die "File $file: $!\n";

my $dbh = DBI->connect("dbi:mysql:db=$arg{D}", $arg{u}, $arg{p});
if (!$dbh) { die "Couldn't connect to $arg{D} using $arg{u} and $arg{p}\n";}

my $q = "SELECT a.feature_id, q.qualifier, a.value"
    . " FROM feature_annotations a, INSDC.qualifier q"
    . " WHERE q.id = a.data_type_id";
my $result = $dbh->selectall_arrayref($q);
my %LU;
my %DATA;
foreach my $row(@$result) {
    my ($id, $qual, $value) = @$row;
    if ($qual eq 'db_xref') {
	my ($k, $v) = split/\:/,$value;
	push @{$DATA{$id}->{$k}}, $v;
	if ($k eq "MGRAST") { $LU{$v} = $id }
    } else {
	push @{$DATA{$id}->{$qual}}, $value;
    }
}    

my $binq = "SELECT m.feature_id, s.name FROM seq_feat_mappings m, sequence_sets s, seq_set_link l WHERE l.seq_id=m.seq_id and s.set_id=l.set_id and s.set_id > 20";
my $bins = $dbh->selectall_hashref($binq, 'feature_id');
foreach my $id (keys %$bins) {
    push @{$DATA{$id}->{'bin'}}, $bins->{$id}->{'name'};
}

while (my $line = <$in>) {
    if ( $line =~ /^>(\S+)/) {
	my $acc = $1;
	my @p;
	foreach my $k (keys %{$DATA{$LU{$acc}}}) {
	    foreach my $v (@{$DATA{$LU{$acc}}->{$k}}) {
		push @p, "$k=\"$v\"";
	    }
	}
	print ">$arg{D}\:$LU{$acc} " . join(" ", @p) . "\n";
    } else {
	$line =~ s/(.{60})/$1\n/g;
	print $line if ($line =~ /\w+/);
    }
}
