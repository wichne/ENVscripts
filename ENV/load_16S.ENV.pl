#!/usr/local/bin/perl

use DBI;
use lib "/home/wcnelson/devel";
use ENV;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use strict;
$| = 1; # Autoflush

######################################## GLOBAL SETTINGS #################

my $usage = " 
    USAGE: $0 -U user -P password -D db [-h -A asmbl_id -V -O prefix -d -T -r]

    REQUIRED:
    -u for username
    -p for password
    -D <database> (e.g. gct)

    OPTIONS:
    -A <asmbl_id> for an old assembly
    -O <RNA prefix> (e.g. CT) default uses first letter genus and species in common..genomes.name
    -T for search results without changes to the database
    -r use previous search results instead of running searches again

    -h print this message and exit
    -V for verbose mode
";

my %args;
&getopts('u:p:d:hV', \%args);
my $VERBOSE = $args{V} ? 1 : 0;
my $server = "orthanc.usc.edu";
my $dbtype = "mysql";
my $dbh = DBI->connect("dbi:$dbtype:server=$server;database=$args{d}", $args{'u'}, $args{'p'});

my $file = $ARGV[0];

######################################## BODY ############################
# DELETE OLD RNA DATA
#foreach my $a (keys %MOLS) { &deleteOldRNAs($dbproc, $a) unless ($DEBUG); } 

open (BTAB, "$file") || die "I can't open '$file': $!\n";
print "\t\tGo through btab output...\n" if ($VERBOSE);
my %INFO;
while (my $match = <BTAB>)  {
    next if ($match =~ /^\#/);
    chomp($match);
    my @info = split(/\t/, $match);
    print "\t\t$info[0] hits $info[4]\n" if ($VERBOSE);
    push @{$INFO{$info[0]}}, \@info;
}

foreach my $id (keys %INFO) {
    my $hits = $INFO{$id};
    for (my $i=0; $i<@$hits; $i++) {
	if (! $hits->[$i]->[0]) { next }
	my ($lo, $hi) = ($hits->[$i]->[2], $hits->[$i]->[3]);
	# determine range that first hit should fall into
	my ($range_lo, $range_hi);
	if ($hits->[$i]->[9] >= 0 ) {
	    $range_lo = $hits->[$i]->[2] - ($hits->[$i]->[7] - 1);
	    $range_hi = $hits->[$i]->[3] + ($hits->[$i]->[6] - $hits->[$i]->[8] + 1);
	} else {
	    $range_lo = $hits->[$i]->[2] - ($hits->[$i]->[6] - $hits->[$i]->[8] + 1);
	    $range_hi = $hits->[$i]->[3] + ($hits->[$i]->[7] - 1);
	}
	if ($range_lo <= 0) { $range_lo = 1}
	my $hit_lo = $hits->[$i]->[7];
	my $hit_hi = $hits->[$i]->[8];
	my ($lo_mod, $hi_mod);
	if ($hits->[$i]->[9] > 0) {
	    $lo_mod = $hits->[$i]->[7] != 1 ? "<" : "";
	    $hi_mod = $hits->[$i]->[8] == $hits->[$i]->[6] ? "" : ">";
	} else {
	    $lo_mod = $hits->[$i]->[8] == $hits->[$i]->[6] ? "" : "<";
	    $hi_mod = $hits->[$i]->[7] != 1 ? ">" : "";
	}

	my @tmp;
	for (my $j=$i+1; $j < @$hits; $j++) {
	    if (&overlaps($hits->[$j]->[2], $hits->[$j]->[3], $range_lo, $range_hi) > 1
		&& $hits->[$i]->[9] eq $hits->[$j]->[9]) {
		$lo = $hits->[$j]->[2] < $lo ? $hits->[$j]->[2] : $lo;
		$hit_lo = $hits->[$j]->[7] < $hit_lo ? $hits->[$j]->[7] : $hit_lo;
		$hi = $hits->[$j]->[3] > $hi ? $hits->[$j]->[3] : $hi;
		$hit_hi = $hits->[$j]->[8] > $hit_hi ? $hits->[$j]->[8] : $hit_hi;
		if ($hits->[$i]->[9] > 0) {
		    $lo_mod = $lo_mod && $hits->[$j]->[7] == 1 ? "" : $lo_mod;
		    $hi_mod = $hi_mod && $hits->[$j]->[8] == $hits->[$j]->[6] ? "" : $hi_mod;
		} else {
		    $lo_mod = $lo_mod && $hits->[$j]->[8] == $hits->[$j]->[6] ? "" : $lo_mod;
		    $hi_mod = $hi_mod && $hits->[$j]->[7] == 1 ? "" : $hi_mod;
		}
		$hits->[$j] = undef;
	    }
	}
	print "$id $lo_mod$lo/$hi_mod$hi to $hit_lo/$hit_hi ($hits->[$i]->[6]) on $hits->[$i]->[9]\n";
	my $featobj = Bio::SeqFeature::Generic->new( 
						     -start        => $lo_mod . $lo, 
						     -end          => $hi_mod . $hi,
						     -strand       => $hits->[$i]->[9], 
						     -primary      => 'rRNA',
						     -display_name => '16S rRNA',
						     );

	my $seq_id = get_seq_id_by_accession($dbh, $id);
	my $feature_id = load_SeqFeature($dbh, $seq_id, $featobj, undef, "0001000");
    }
}


sub overlaps {
    my ($lo, $mid1, $mid2, $hi) = sort {$a<=>$b} @_;
    return ((abs($_[1] - $_[0]) + 1) + (abs($_[3] - $_[2]) + 1)) - ($hi - $lo + 1);
}
