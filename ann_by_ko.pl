#!/usr/bin/perl

#Edited by Jmo 11/09/2015 to use ghostKO/Koala KEGG files
#Edited by Jmo 07/18/14 to change in and parse the feature_acc to accept a KAAS file with multiple acc 

use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use strict;

my %arg;
&getopts('D:i:u:p:', \%arg);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $user = $arg{u} ? $arg{u} : $ENV{USER};
my $password = $arg{p} or die "No password provided.";
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg{D}", $user, $password);

my $ko_file = $arg{i};

my $ko_q = "SELECT * FROM egad.ko";
my $koref = $dbh->selectall_hashref($ko_q, 'ko');

if ($ko_file) {
    print STDERR "Loading data from $ko_file...\n";
    open my $in, $ko_file or die "Can't open $ko_file: $!\n";
    
    ## 1. Read in line
    while (my $line = <$in>) {
	next if ($line =~ /^#/);
	chomp $line;
	my ($feat_acc,
	    $ko) = split/\s+/, $line;
	&load_ko_annotation($dbh, &get_feature_id_by_accession($dbh,$feat_acc), $ko) if ($ko);
    }

} else {
    # Get ko hit info from database
    print STDERR "Retrieving KO hit data from database $arg{D}...\n";
    
    my $feat_count_q = "SELECT count(feature_id) from sequence_features where feat_type='CDS'";
    my $ko_feat_count_q = "SELECT count(distinct feature_id) from feature_evidence where ev_type='KO'";
    my ($fc) = $dbh->selectrow_array($feat_count_q);
    my ($hfc) = $dbh->selectrow_array($ko_feat_count_q);
    print STDERR "I see $hfc features with ko hits out of $fc features total.\nDoes that sound right?(Y/n) ";
    my $answer = <STDIN>;
    if ($answer =~ /^n/i) { print STDERR " Dying! "; die; }

    my $ko_q = "SELECT feature_id, ev_accession, score"
	. " FROM feature_evidence"
	. " WHERE ev_type='KO'";
    print "$ko_q\n";
    my $ko_ev_ref = $dbh->selectall_arrayref($ko_q);
    foreach my $row (@$ko_ev_ref) {
	my $ko_acc = $row->[1];
	print "$ko_acc\n";
	&load_ko_annotation($dbh, $row->[0], $ko_acc);
    }
}


sub load_ko_annotation {
    my $dbh = shift;
    my $feat_id = shift;
    my $ko_acc = shift;

    # clear out existing data
    my $del_q = "DELETE FROM feature_annotations"
	. " WHERE feature_id=$feat_id"
	. " AND rank = 5"
	. " AND source = 'ghostKOALA'";
    $dbh->do($del_q);
    
    my $upd_q = "INSERT feature_annotations"
	. " (feature_id, data_type_id, value, rank, source, edit_by)"
	. " VALUES ($feat_id, ?, ?, 5,\"ghostKOALA\", USER())"
	;
    my $sth = $dbh->prepare($upd_q);
    if ($koref->{$ko_acc}->{product}) {
	$sth->execute(66, $koref->{$ko_acc}->{product});
    }
    if ($koref->{$ko_acc}->{gene_sym}) {
	$sth->execute(35, $koref->{$ko_acc}->{gene_sym});
	}
    if ($koref->{$ko_acc}->{ec}) {
	$sth->execute(1, $koref->{$ko_acc}->{ec});
    }
}

exit;
