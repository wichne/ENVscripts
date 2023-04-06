#!/usr/bin/perl

# This program deletes a feature from all tables

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;
$| = 1;

my %arg;
&getopts('D:u:p:f:l:', \%arg);
my $dbh = &connect(\%arg);

my @ids;

if ($arg{'l'}) {
    if (-e $arg{'l'}) {
        open(my $list, $arg{'l'}) or die "Can't open $arg{l} for read: $!\n";
        while (my $fid = <$list>) {
            chomp $fid;
            push @ids, $fid;
        }
    }
} elsif ($arg{'f'}) {
    push @ids, $arg{'f'};
}

delete_sequence_features($dbh, @ids);
