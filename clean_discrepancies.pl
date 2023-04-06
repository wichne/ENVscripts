#!/usr/bin/perl

# This program tries to automagically fix a number of the 
# common complaints in a tbl2asn discrepancy file including
# EC in name, TC in name, common bad characters and strings
# more

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use Getopt::Std;
$| = 1;

my %arg;
&getopts('D:u:p:n:i:', \%arg);
my $db = $arg{'D'};
my $dbh = &connect(\%arg);

#system("move_EC_from_descriptor.pl -D $db -u $arg{'u'} -p $arg{'p'}");
#system("move_TC_from_descriptor.pl -D $db -u $arg{'u'} -p $arg{'p'}");

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
    ;

my $result = $dbh->selectall_arrayref($q);

my $val_u = "UPDATE feature_annotations set value = ?"
    . " WHERE feature_id=?"
    . " AND value = ?";
my $val_sth = $dbh->prepare($val_u);

my $ec_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, ann_rank, source)"
    . " VALUES (?, 1, ?, ?, ?)";
my $ec_sth = $dbh->prepare($ec_i);

my $tc_i = "INSERT INTO feature_annotations (feature_id, data_type_id, value, ann_rank, source)"
    . " VALUES (?, 98, ?, ?, ?)";
my $tc_sth = $dbh->prepare($tc_i);

foreach my $row(@$result) {
    my ($fid, $val, $rank, $source) = @$row;
 
    # [EC:3.1.3.97] or (EC 2.1.1.61) or [EC:2.1.1.61 1.5.-.-]
    print "$fid :: $val ->\n";
    my $new_val = $val;
    if ($new_val =~ /[\(\[]\s*EC\:?\s*([^\)\]]+)\s*[\)\]]/) {
        my @ecs;
        while ($new_val =~ /[\(\[]\s*EC\:?\s*([^\)\]]+)\s*[\)\]]/g) {
            my $whole_string = $&;
            push @ecs, split(/[, ]+/, $1);
            $new_val =~ s/\Q$whole_string\E//;
            #$new_val =~ s/\(\)//;
        }
        print "\t$new_val\n";

        foreach my $ec (@ecs) {
	        if ($ec =~ /\d+(\.(\d+|\-)){1,2}(\.[nB]?\d+)?/) {
	            print "\t-> $ec\n";
                $ec_sth->execute($fid, $ec, $rank, $source);
	        } else {
	            print "bad ec: $ec\n";
	            next;
	        }
        }
    }

    # Handle TC in value
    if ($new_val =~ /\s*\(\s*TC \d+\.[A-Z](\.\d+){3}\s*\)/) {
        my @tcs;
        while ($new_val =~ /\s*\(\s*TC (\d+\.[A-Z](\.\d+){3})\s*\)/g) {
            my $whole_string = $&;
            push @tcs, split(/[, ]+/, $1);
            $new_val =~ s/\Q$whole_string\E//;
            #$new_val =~ s/\(\)//;
        }
        print "\t$new_val\n";

        foreach my $tc (@tcs) {
	        if ($tc =~ /\d+\.[A-Z](\.\d+){3}/) {
	            print "\t-> $tc\n";
	            $tc_sth->execute($fid, $tc, $rank, $source);
        	} else {
	            print "bad tc: $tc\n";
	            next;
	        }
        }
   }

# if name like " and related (.*(enzyme|protein|\w*ase))s" switch to "-like protein"
# or if there is the descriptor change it to the captured descriptor
    if ($new_val =~ /( and related (.*(enzyme|protein|\w+ase))s)/) {
        #print "$new_val\t->\t";
        #if ($3 ne $2) {
        #    $new_val = "${3}-related protein";
        #} else {
            $new_val =~ s/\Q$1\E/-related protein/;
        #}
        print "\t$new_val\n";
    }

# if name like " \w*ases\b | enzymes | proteins" make singular
    if ($new_val =~  /((\w+ase|enzyme|protein)s)\b/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E/\Q$2\E/;
        print "\t$new_val\n";
    }

# if name like " > .*" delete from > on.
   if ($new_val =~  / > .*/) {
        #print "$new_val\t->\t";
        $new_val =~ s/ > .*//;
        print "\t$new_val\n";
    }

# if name like "hypothetical protein.+" change to hypothetical protein
    if ($new_val =~  /^(hypothetical protein)(.+)/i) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$2\E//;
        print "\t$new_val\n";
    }

# if name like "fragment( of )?"i
    if ($new_val =~  /(\s*fragment( of )?)/i) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E//;
        print "\t$new_val\n";
    }

# if name like "FIG\d{8}(: )?" delete
    if ($new_val =~  /(\s*(FIG(fam)?\d{6,8}|FOG)(: )?)/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E//;
        print "\t$new_val\n";
    }

# if name like "\w+\_\w*\d+" delete
    if ($new_val =~  /(\w+\_\w*\d+)\b/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E//;
        print "$new_val\n";
    }
# if name like "\?" replace with ", putative"
    if ($new_val =~  /\?/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\?/, putative/;
        print "\t$new_val\n";
    }
# if name like "\@" replace with "or"?
    while ($new_val =~  /\s\@\s/g) {
        #print "$new_val\t->\t";
        $new_val =~ s/\@/or/;
        print "\t$new_val\n";
    }
# if name like "Chromosome undetermined.*" replace with "hypothetical protein"
    if ($new_val =~  /(Chromosome undetermined.*)/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E/hypothetical protein/;
        print "\t$new_val\n";
    }
# if name like " homolog" replace with "-related"
    if ($new_val =~  /( homolog)/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E/-related/;
        print "\t$new_val\n";
    }

    if ($new_val =~  /(\s*involved in .*)/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E//;
        print "\t$new_val\n";
    }

    # replace ';' with '/'
    if ($new_val =~  /(\s*\;\s*)/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E/ \/ /g;
        print "\t$new_val\n";
    }

    if ($new_val =~  /([\_\,\/]\s*)$/) {
        #print "$new_val\t->\t";
        $new_val =~ s/${1}$//;
        print "\t$new_val\n";
    }

    if ($new_val =~  /\b(likely|possible|possibly|probably)\b/) {
        #print "$new_val\t->\t";
        $new_val =~ s/\Q$1\E/putative/;
        print "\t$new_val\n";
    }

    # Get rid of protein accessions
    if ($new_val =~ /[A-Z][A-Za-z]{1,2}\d{4,5}/) {
        while ($new_val =~ /([A-Z][A-Za-z]{1,2}\d{4,5})/g) {
            my $acc = $1;
            unless ($acc =~ /^(DUF|UPF|COG)/) {
                $new_val =~ s/\s*$acc//;
                print "\t$new_val\n";
            }
        }
    }

    if ($new_val =~ /(\(\s*\))/) {
        $new_val =~ s/\Q$1\E//;
        print "\t$new_val\n";
    }

    $val_sth->execute($new_val, $fid, $val) if ($val ne $new_val);
} 