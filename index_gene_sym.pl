#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

my %args;
&getopts('D:p:u:n:i:', \%args);

if (!($args{D} && $args{p} && ($args{n} || $args{i}))) {
    die "USAGE:
index_gene_sym.pl -D db -p dbpswd [-u user] [-n set_name or -i set_id]
";
}

my $dbh = connect(\%args);

my $set_id = $args{i} ? $args{i} : set_name_to_id($dbh, $args{n});
my $feat_ref = get_seq_features_by_set_id($dbh, $set_id);
foreach my $featobj(keys %$feat_ref) {
}


my $the_query = "SELECT a.feature_id, a.rank, a.source, a.value,"
    . " s.seq_id, s.seq_length, m.feat_min"
    . " FROM feature_annotations a, seq_feat_mappings m, seq_set_link l, sequences s"
    . " WHERE a.data_type_id=35"
    . " AND m.feature_id=a.feature_id"
    . " AND l.seq_id = m.seq_id"
    . " AND s.seq_id=m.seq_id"
    . " AND l.set_id=$set_id"
    . " ORDER BY a.feature_id, a.rank";

my $res = $dbh->selectall_arrayref($the_query);

my $last;
foreach my $result(@$res) {
    my ($fid, $rank, $source, $value, $sid, $len, $feat_min) = @$result;
    if ($fid == $last)  {next}
    push @{$DATA{$value}}, $result;
    $last = $fid;
}

foreach my $gene (keys %DATA) {
    if (@{$DATA{$gene}} > 1) {
	my $n = 1;
	foreach my $fobj (sort {$b->[5] <=> $a->[5] ||
				    $a->[4] <=> $b->[4] ||
				    $a->[6] <=> $b->[6]} @{$DATA{$gene}}) {
	    my ($fid, $rank, $source, $value, $sid, $len, $feat_min) = @$fobj;
	    my $new_value = $value . "-" . $n++;
	    print "$fid\t$rank\t$source\t$new_value\n";
	    my $update_q = "UPDATE feature_annotations set value=\"$new_value\""
		. " WHERE feature_id=$fid AND rank = $rank AND source = \"$source\""
		. " AND data_type_id=35 AND value=\"$value\"";
	    $dbh->do($update_q);
	}
    }
}
