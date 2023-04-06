#!/usr/bin/perl

use Getopt::Std;
use DBI;
use strict;
use lib $ENV{SCRIPTS};
use ENV;

my %arg;
&getopts('D:u:p:i:', \%arg);

my $user = $arg{'u'};
my $password = $arg{'p'};
my $db = $arg{'D'};
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}

my $infile = $arg{'i'};
if (! $infile) { die "No prodigal file provided with -i\n";}
if (! -r $infile) { die "$infile is not readable: $!\n";}

open my $in, $infile or die "Can't open $infile: $!\n";
## prodigal output:
# DEFINITION  seqnum=1;seqlen=1073940;seqhdr="contig-100_0 length_1073940 read_count_5843341";version=Prodigal.v2.60;run_type=Metagenomic;model="15|Chlorobium_phaeobacteroides_BS1|B|48.9|11|0";gc_cont=48.90;transl_table=11;uses_sd=0
#FEATURES             Location/Qualifiers
#     CDS             <2..772
#                     /note="ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.554;conf=100.00;score=81.72;cscore=78.50;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;"
# etc...
# //
 
while (my $line = <$in>) {
    chomp $line;
    if ($line =~ ^/\/\//) { &process_record($record) if ($record) }
    else { $record
}
