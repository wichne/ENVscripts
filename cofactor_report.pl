#!/usr/bin/perl
use DBI;
use lib "/share/scripts";
use ENV;
use Getopt::Std;
use strict;

my %arg;
&getopts('D:u:p:',\%arg);
my $db = $arg{D};
my $dbp= $arg{p};
my $dbu = $arg{u} ? $arg{u} : $ENV{USER};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $dbu, $dbp);

my $MR = &get_mainroles($dbh);
my $SR = &get_subroles($dbh);

my $setref = get_seq_sets($dbh);
foreach my $setid (sort {$a<=>$b} keys %$setref) {
    my $featref = get_seq_features_by_set_id($dbh, $setid, "CDS");
    foreach my $featid (sort {$a<=>$b} keys %$featref) {
	my ($acc, $role, $product, $ec, $hmmcof, $eccof);
	foreach my $a (@{$featref->{$featid}->{'accessions'}}) {
	    if (!$acc) { $acc = $a }
	    elsif ($a =~ /^ITZ/) { $acc = $a }
	}
	$hmmcof = join(" ", @{$featref->{$featid}->{'cofactor'}->{'hmm'}}) if (defined $featref->{$featid}->{'cofactor'}->{'hmm'});
	$eccof = join(" ", @{$featref->{$featid}->{'cofactor'}->{'ec'}}) if (defined $featref->{$featid}->{'cofactor'}->{'ec'});
	
	if (defined $featref->{$featid}->{role}->{mainrole}) {
	    for (my $i=0; $i<@{$featref->{$featid}->{role}->{mainrole}}; $i++) {
		$role .= "$MR->{$featref->{$featid}->{role}->{mainrole}->[$i]}->{name}::$SR->{$featref->{$featid}->{role}->{subrole}->[$i]}->{name};";
	    }
	}
	my %ANNOT;
	foreach my $key ('product', 'gene', 'EC_number') {
	    my ($best_rank) = sort {$a<=>$b} keys %{$featref->{$featid}->{'annotation'}->{$key}};
	    $ANNOT{$key} = join(" ", sort (map $_->{'value'}, @{$featref->{$featid}->{'annotation'}->{$key}->{$best_rank}}));
	}
#	my ($best_rank) = sort {$a<=>$b} keys %{$featref->{$featid}->{'annotation'}->{product}};
#	$product = $featref->{$featid}->{'annotation'}->{product}->{$best_rank}->{'value'};
	
	printf "$featid\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", ($setref->{$setid}->{'name'}, $acc, $ANNOT{'product'}, $ANNOT{'gene'}, $ANNOT{'EC_number'},$role, $hmmcof, $eccof);
    }
}
exit;
