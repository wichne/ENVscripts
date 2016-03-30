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

my %DATA;
my %COF;
my $cofactor_q = "SELECT s.name, cofactor, count(distinct a.feature_id)"
    . " FROM feature_annotations a, seq_feat_mappings m, sequence_sets s, seq_set_link l, egad.ec_cofactor_link e"
    . " WHERE s.set_id=l.set_id and m.seq_id=l.seq_id"
    . " AND a.feature_id=m.feature_id and a.data_type_id=1"
    . " AND a.value=e.ec"
    . " GROUP by s.name, cofactor";
my $sth = $dbh->prepare($cofactor_q);
$sth->execute();
while (my @row = $sth->fetchrow_array) {
    $DATA{$row[0]}->{$row[1]} = $row[2];
    $COF{$row[1]} = 1;
}

print "Genome\t";
print join("\t", sort keys %COF),"\n";

foreach my $genome (sort keys %DATA) {
    print "$genome";
    foreach my $cof (sort keys %COF) {
	print "\t$DATA{$genome}->{$cof}";
    }
    print "\n";
}

my %VHMMs = {'Thiamine' => ['PF01558', 'PF01855', 'TIGR03710'],
	     'TPP' => ['PF00456', 'PF00676', 'PF02775', 'PF02776', 'PF02775', 'PF13292', 'TIGR00118', 'TIGR00204', 'TIGR00232', 'TIGR00239', 'TIGR02176', 'TIGR03182'],
	     'FAD' => ['PF00070','PF00890','PF01266','PF01494','PF01565','PF03441','PF04898','PF07992'],
	     'FMN' => ['PF01070','PF01180','PF01243','PF01645','PF12766','TIGR01036'],
	     'flavin' => ['PF00070','PF00146', 'PF00329','PF00346','PF00361','PF00420','PF00499','PF00507','PF00662','PF00970','PF01058','PF01264','PF01512','PF01593','PF02219','TIGR00033','TIGR00179','TIGR00562','TIGR01350','TIGR01551','TIGR01770','TIGR01811','TIGR01957','TIGR01961','TIGR01962','TIGR01972','TIGR01974','TIGR02151'],
	     'NAD' => ['PF00175','PF00465','PF00670','PF01370','PF01761','PF02826','PF04321','PF05221','PF13450','PF13602','PF13685','TIGR00936','TIGR01181'],
	     'NAD(P)' => ['PF00106','PF01073','PF02882'],
	     'NADH' => ['PF00070','PF00107','PF07992'],
	     'NADPH' => ['PF08659'],

exit;
