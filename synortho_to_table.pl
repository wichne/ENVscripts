#!/usr/bin/perl
#This is a script that makes output from the 
use strict;
use DBI;
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;
$| = 1;

my %arg;
&getopts('D:u:p:i:o:',\%arg);
my $password = $arg{p} or die "Need to provide password with -p\n";
my $db = $arg{D};
my $user = $arg{u} ? $arg{u} : $ENV{USER};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db",
		       $user, $password, {'AutoCommit'=> 0});
# open file
my $infile = $arg{i} or die "Please provide input file with -i";
open my $in, $infile or die "Can't open '$infile': $!\n";

my %SETS;
my %D;
while (my $line = <$in>) {
    chomp $line;
    my ($id, @acc) = split/\t/, $line;
    my $synid = sprintf "SYNORTHO%05d", $id;
    foreach my $a(@acc) {
	$a =~ /([^\_]+)\_/;
	my $set = $1;
	push @{$D{$synid}{$set}}, $a;
	$SETS{$set}++;
    }
}

print "synorthoid\t" . join("\t", sort keys %SETS) . "\n";
    
foreach my $synid (sort keys %D) {
    my ($product, $ec, $mainrole);
    print "$synid";
    foreach my $set (sort keys %SETS) {
	if (defined $D{$synid}{$set}) {
	    my @acc = @{$D{$synid}{$set}};
	    my ($fid) = get_feature_id_by_accession($dbh,$acc[0]);
	    my $featref = get_features_by_feature_id($dbh, $fid);
	    if (! $product && defined $featref->{$fid}->{'annotation'}) {
		my ($k) = sort {$a<=>$b} keys %{$featref->{$fid}->{'annotation'}};
		my ($source, $ref) = each %{$featref->{$fid}->{'annotation'}->{$k}};
		$product = $ref->{'product'}->[0];
	    }
	    if (! $ec && defined $featref->{$fid}->{'annotation'}) {
		my ($k) = sort {$a<=>$b} keys %{$featref->{$fid}->{'annotation'}};
		my ($source, $ref) = each %{$featref->{$fid}->{'annotation'}->{$k}};
		$ec = $ref->{'ec'}->[0];
	    }

	    if (! $mainrole) {
		my $roleref = get_ev_role_by_feature_id($dbh, $fid);
		$mainrole = join(" ", @{$roleref->{'mainrole'}});
	    }
	    print "\t\"" . join(" ", @acc) . "\"";
	} else {
	    print "\t";
	}
    }
    print "\t$product\t$ec\t$mainrole\n";
}
