#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;

my %arg;
&getopts('D:u:p:', \%arg);
my $db = $arg{'D'};
my $dbh = &connect(\%arg);

my $q = "SELECT feature_id, value, ann_rank, source"
    . " FROM feature_annotations"
    . " WHERE data_type_id=66"
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
