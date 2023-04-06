#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;

my %arg;
&getopts('D:u:p:i:n:', \%arg);
my $db = $arg{'D'};
my $dbh = &connect(\%arg);
my $set_id;
if ($arg{'n'}) {
    $set_id = &set_name_to_id($dbh, $arg{'n'});
} elsif ($arg{'i'}) {
    $set_id = $arg{'i'};
} else { die "Barf;\n"; }

my $q = "SELECT a.feature_id, a.value, a.ann_rank, a.source"
    . " FROM feature_annotations a, seq_feat_mappings m, seq_set_link l"
    . " WHERE data_type_id=66"
    . " AND l.set_id = $set_id"
    . " AND m.seq_id=l.seq_id"
    . " AND a.feature_id=m.feature_id"
    . " AND value REGEXP '[\\(\\\\[] ?EC[ :]'"
#    . " AND value REGEXP '[\\(\\\\[]EC[ :][0-9]+(\.[0-9\-]+){3}[\\)\\\\]]'"
#    . " AND value REGEXP '[\\\\[]EC[ :][0-9]+(\.[0-9\-]+){3}[\\\\]]'"
    ;

my $result = $dbh->selectall_arrayref($q);

my $val_u = "UPDATE feature_annotations set value = ?"
    . " WHERE feature_id=?"
    . " AND value = ?";
my $val_sth = $dbh->prepare($val_u);
my $ec_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, ann_rank, source)"
    . " VALUES (?, 1, ?, ?, ?)";
my $ec_sth = $dbh->prepare($ec_i);

foreach my $row(@$result) {
    my ($fid, $val, $rank, $source) = @$row;
    my @ecs;
    my $new_val = $val;
    while ($val =~ /[\(\[] ?EC[ \:]([^\)\]]+)[\)\]]/g) {
        my $whole_string = $&;
        push @ecs, split(/[, ]+/, $1);
        $new_val =~ s/\Q$whole_string\E//;
        #$new_val =~ s/\(\)//;
    }
    print "$fid :: $val ->\n";
    print "\t$new_val\n";
    $val_sth->execute($new_val, $fid, $val);

    foreach my $ec (@ecs) {
	if ($ec =~ /\d+(\.(\d+|\-)){1,2}(\.[nB]?\d+)?/) {
	    print "\t$ec\n";
	    $ec_sth->execute($fid, $ec, $rank, $source);
	} else {
	    print "bad ec: $ec\n";
	    next;
	}
    }
}
