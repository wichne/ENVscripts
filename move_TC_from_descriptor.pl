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
    . " AND value REGEXP '\\(TC [0-9]+\.[A-Z](\.[0-9]+){3}\\)'"
#    . " AND value REGEXP '[\\(\\\\[]EC[ :][0-9]+(\.[0-9\-]+){3}[\\)\\\\]]'"
#    . " AND value REGEXP '[\\\\[]EC[ :][0-9]+(\.[0-9\-]+){3}[\\\\]]'"
    ;

my $result = $dbh->selectall_arrayref($q);

my $val_u = "UPDATE feature_annotations set value = ?"
    . " WHERE feature_id=?"
    . " AND value = ?";
my $val_sth = $dbh->prepare($val_u);
my $tc_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, ann_rank, source)"
    . " VALUES (?, 98, ?, ?, ?)";
my $tc_sth = $dbh->prepare($tc_i);

foreach my $row(@$result) {
    my ($fid, $val, $rank, $source) = @$row;
    my @tcs;
    my $new_val = $val;
    while ($val =~ /[\(\[] ?TC[ \:]([^\)\]]+)[\)\]]/g) {
        my $whole_string = $&;
        push @tcs, split(/[, ]+/, $1);
        $new_val =~ s/\Q$whole_string\E//;
        #$new_val =~ s/\(\)//;
    }
    print "$fid :: $val ->\n";
    print "\t$new_val\n";
    $val_sth->execute($new_val, $fid, $val);

    foreach my $tc (@tcs) {
	if ($tc =~ /\d+\.[A-Z](\.\d+){3}/) {
	    print "\t$tc\n";
	    $tc_sth->execute($fid, $tc, $rank, $source);
	} else {
	    print "bad tc: $tc\n";
	    next;
	}
    }
}
