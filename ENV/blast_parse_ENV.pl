#!/usr/bin/perl

use Cwd;
use DBI;
use Getopt::Std;
use TIGR::FASTArecord;
use TIGR::FASTAreader;
use Bio::Annotation::Collection;
use strict;
#use lib "/home/wcnelson/devel/";
use vars qw(%ber_rank);
require "autoAnnotate.data";
use ENV;

$| = 1;
my $VERBOSE=0;
my $DEBUG;

print STDERR localtime(time) . "\n";
my @filtertext;
my %opt;
&getopts('P:D:f:U:', \%opt);

my $dbh = DBI->connect("dbi:mysql:server=localhost;database=$opt{D}", $opt{U}, $opt{P});

my $per_id_cutoff = 30;
my $length_cutoff = 80;
my $annotation_ref = {};
&parse_BLAST_tab($dbh, $annotation_ref, $opt{f});

while (my ($feature_id, $ref) = each (%$annotation_ref)) {
    my $annObj = Bio::Annotation::Collection->new();
    foreach my $key("product", "gene", "EC_number") {
	if (defined $ref->{$key}) {
	    my $sv = Bio::Annotation::SimpleValue->new(-value=>$ref->{$key});
	    $annObj->add_Annotation($key, $sv);
	}
    }
    &delete_feature_annotations($dbh, $feature_id);
    &load_feature_annotations($dbh, $feature_id, $annObj);
}

exit();

sub parse_BLAST_tab {
    my($dbproc, $ref, $file) = @_;

    ## tab records considered:
    ## [0] query accession
    ## [1] length of query
    ## [2] query start
    ## [3] query end
    ## [4] match accession
    ## [5] match description
    ## [6] match length
    ## [7] match hit start 
    ## [8] match hit end
    ## [9] strand (1/-1) or frame (1/2/3/-1/-2/-3)
    ## [10] bit score
    ## [11] %similarity
    ## [12] %identity
    ## [13] E-value

    if (-r $file) { # open btab only if it is readable
	my $TAB = {};
	my $best_pvalue;
	my $no_count;
	my $last_acc;
	open(TAB, $file);
	while (my $line = <TAB>) { # cycle thru btab
	    if ($line =~ /^\#/) { next }
	    chomp($line);
	    my $tab_ref = {};
	    my @tab =  split(/\t/,$line); # make array of btab info for the row
	    if ($tab[0] ne $last_acc) {
		foreach my $x (sort { $TAB->{$a}->[13] <=> $TAB->{$b}->[13] }keys %$TAB) {
		    my $tab = $TAB->{$x};
		    $$tab[5] =~ s/\s+/ /g; # replace multiple whitespace with a single space
		    
		    # here we want to compare the current line with the previous best line.
		    # TIGR > Characterized > OMNI > EGAD > Anything else
		    # So replace if info exists and is better.
		    
		    # Get best annotation from header
		    my $tab_ref = &parse_header($dbproc, $$tab[5]);
		    $tab_ref->{product} = &add_annotation_level($tab_ref->{product}, $$tab[12], $$tab[1]);
		    
		    # Determine hit rank
		    $tab_ref->{acc} =~ /^(.*?)\|(\w+)/; # everything before the first | which is the db
		    my $source = $1;
		    my $acc = $2;
		    if ($source eq "OMNI" && $acc !~ /^NT/) { $source = "TIGR" }
		    elsif ($tab_ref->{exp} == 1) { $source = "characterized" }
		    $tab_ref->{rank} = defined $ber_rank{$source} ? $ber_rank{$source} : 100;
		    print "\tdb source: $source: rank $tab_ref->{rank}\n" if ($DEBUG);
		    # This whole ranking thing also has to be subject to some hit-length check
		    if ($$tab[12] < $per_id_cutoff) { $tab_ref->{rank} += 6 }
		    $tab_ref->{length} = $$tab[14];
		    $tab_ref->{perid} = $$tab[12];
		    if ($tab_ref->{length} * 100 <= $length_cutoff) { $tab_ref->{rank} += 6 }
		    
		    # If this hit ranks, compare its annotation to the current best set
		    if ($tab_ref->{rank} <= $ref->{$$tab[0]}->{rank} ||
			! defined $ref->{$$tab[0]}->{rank}) {
			print "!BER hit $tab_ref->{rank} meets or beats current $ref->{rank}\n" if ($VERBOSE);
			# this sub compares ident info, and sets best info.
			&compare_info($ref->{$$tab[0]}, $tab_ref);
		    }
		}
		$no_count = 0; # reset no_count so it only counts consecutive no-hits.
		$best_pvalue = undef;
		$TAB = {};
	    }

	    if ($tab[13] == 0) { $tab[13] = 1e-320 }
	    if (! defined $best_pvalue) { $best_pvalue = $tab[13] }

	    print "--------------------------\n$tab[4]\n" if ($VERBOSE);
#	    if (&filtercheck($tab[5])) { next };

	    # this controls for genes at end of molecule that don't get full extend
	    my $gene_length = $tab[1];
	    my $tab_length = ($tab[3] - $tab[2] + 1) / $gene_length; # fraction length
	    print "LENGTH: " . $tab_length * 100 . "%\n" if ($DEBUG);

	    # Check some score stuff.
	    if ($tab[12] > 30 &&
		# believe it or not, there is empirical data behind this cutoff
		$tab[10] * $tab_length > 4 &&
		# this tries to see how close the current best score is to this score
		# looks at the log of the negative log of the pvalue
		# ie a pvalue of 1e-100 would evaluate to 2, so take everything down to 1e-1e1.5 which is 2.38e-32.
		# there is no math or science behind this cutoff - for now.
		($best_pvalue && $tab[13] &&
		 log(-log($best_pvalue)/log(10))/log(10) - log(-log($tab[13])/log(10))/log(10) < 0.5) ) {
		print "\tMeets score requirements: ($tab[10]) ($tab[12], $tab_length) ($tab[13])\n" if ($VERBOSE);
		$TAB->{$tab[4]} = \@tab;
		$TAB->{$tab[4]}->[5] = $tab[4] . " " . $tab[5];
		$TAB->{$tab[4]}->[14] = $tab_length;
	    } else {
		$no_count++;
		print "\tNO COUNT: $no_count\n" if ($VERBOSE);
#		if ($no_count > 5) {last;} # breaks out of while?
	    }
	    # Don't want to set best pvalue to identity hit since that might be hard to match.
	    $best_pvalue = undef if ($tab[12] == 100 && $tab_length >= 99);
	    $last_acc = $tab[0];
	}
	close(TAB);
    } 
    else {
	print LOG "Could not open/read $file\n";
    } ## if (-e $file)
}

sub filtercheck {
    my ($string) = @_;
    my $filter = 0;
    # @filtertext should be global variable
    foreach my $word(@filtertext) {
	if ($string =~ /$word/) { $filter = 1}
    }
    return($filter);
}
    
sub parse_header {
    my ($dbproc, $header) = @_;
    my (@x, @y, @z, $name, $species, $gene, $temp, $query, $x, $p_id, @ec, $genus, $exp);
    my $best = {};
    print "parsing header: $header\n" if ($DEBUG);

    $header =~ s/\s+/ /g;
#    @x = split(/\s*\^\|\^/,$header);  # the seperator between panda records.
    @x = split(/\s\>/,$header);  # the seperator between genbank records.
    foreach my $nr_entry (@x) {
	# split accesion and description
	my $test = {};
	$nr_entry =~ /^(\S+)\s+(.*)/;
	$test->{acc} = $1;
	my $desc = $2;

	print "INDIV ACC : $desc\n" if ($DEBUG);

	# grab panda metadata
	if ($desc =~ s/ \(exp(erimental)?=(\-?\d+); wgp=(\-?\d+); cg=(\-?\d+); \)//) {
	    $test->{exp} = $2;
	    $test->{wgp} = $3;
	    $test->{cg} = $4;
	}

	# Also peek in our characterized table.
	if (!$test->{exp}) { $test->{exp} = &checkForExperimental($dbproc, $test->{acc}) }

	# if we already have an experimentally characterized match, don't bother with the ones that aren't.
	if ($best->{exp} == 1 && $test->{exp} != 1) { next }

	# find species
#	if ($desc =~ s/\{([^;\}]+);.*\}//) { $test->{species} = $1 }
	if ($desc =~ s/\[([^\]]+)\]//) { $test->{species} = $1 }
	print "found species: $test->{species}\n" if ($VERBOSE);

	# Maybe what I should actually be doing here is using the Swissknife (SWISS::Entry)
	# module (Sourceforge) to parse info from our local flatfile.
	# Easier, more accurate. At least we could get gene out.
	if ($test->{acc} =~ /^SP/) {
	    while ($desc =~ s/\(EC ([^\)]+)\)//) { 
		push @ec, $1;
	    }
	    $test->{EC_number} = join " ", @ec;
	    # SwissProt has alternate names (synonyms in SP talk) in parens following the EC
	    while ($desc =~ /(\s*\([^\)]*\)\s*\.?\s*)$/) { 
		$desc =~ s/\s*\([^\)]*\)\s*\.?\s*$//;
	    }
	    # bifunctional proteins have separate functions enclosed in brackets. kill them.
	    $desc =~ s/\[(includes|contains)\:.*\]//i;
	    $test->{product} = $desc;
	}
	
	elsif ($test->{acc} =~ /^PIR/) {
	    # this may not even exist anymore with the advent of Uniprot
	    next;
	    $desc =~ /^(.*) - .*$/;
	    $test->{product} = $1;
	    if ( $test->{product} =~ s/\(EC ([^\)]+)\)//) {
		$test->{EC_number} = $1;
	    } elsif ( $test->{product} =~ s/\(([\w\-\/\d]{3,})\)// ) {
		$test->{gene} = $1;
	    }
	    $test->{product} =~ s/ ?\[.*//;
	} elsif ($test->{acc} =~/OMNI/) {
 	    my($x, $p_id) = split(/\|/, $test->{acc});
 	    # why are we doing this instead of taking it from the btab? panda header doesn't seem to have complete info
	    # because we get more info more accurately
	    my $dbh = DBI->connect("dbi:Sybase:server=$ENV{DSQUERY}", "access", "access");
 	    my $table = ($p_id =~ /^NT/) ? "omnium..nt_ident" : "omnium..ident";
 	    $query = "select i.com_name, i.gene_sym, i.ec_num, t.genus, t.species"
 		. " from $table i, omnium..asm_feature a, omnium..asmbl_data d, omnium..taxon_link tl, omnium..taxon t "
 		. " where i.locus = \"$p_id\""
 		. " and i.locus=a.locus"
		. " and a.asmbl_id=d.asmbl_id"
		. " and d.db_data_id=tl.db_taxonl_id"
		. " and tl.taxon_uid=t.uid";
 	    print "$query\n" if($DEBUG);
 	    ($test->{product}, $test->{gene}, $test->{EC_number}, $test->{genus}, $test->{species}) = $dbh->selectrow_array($query);
	    $dbh->disconnect;
 	    ## species not currently used
	    $test->{species} = "$test->{genus} $test->{species}";
	}

	else { # ie if GB
	    # FRAMESHIFT and the like sometimes appear between the gene
	    # symbol and the species. This needs to be removed
	    $desc =~ s/(authentic )?frameshift//i;
	    $desc =~ s/(authentic )?point mutation//i;
	    my $gene;
	    if ($desc =~ /\(([a-z]{2,4}[A-Z\d\-]{,3})\)/) {
		$gene = $1;
		$desc =~ s/\($gene\)\s*//;
	    }

	    # get rid of trailing stuff in parens
	    $desc =~ s/\(.*\)\s*$//;
	    
	    ($test->{gene}, $test->{product}) = ($gene, $desc);
	}

	&clean_annotation($test);
	&compare_info($best, $test); 
    }
    return ($best);
}

sub clean_annotation {
    my $ref = shift;
    # this subroutine should hold ALL error checking for annotation info
    my ($name, $gene, $ec) = (\$ref->{product}, \$ref->{gene}, \$ref->{EC_number});

    print "parse_BER_name: $$name / $$gene / $$ec\n" if ($DEBUG);

    print LOG "Cleaning '$$name'.....\n";

    # this is to try to get rid of accessions that may be part of the name
    $$name =~ s/\b[A-Z]{2,}_?[A-Z]*[\d\.]{3,}[a-z]?\s*//g;
    $$name =~ s/\b[A-Z][a-z]{1,2}\d*\-?\d{3,}[a-z]?//g;
    
    # try to identify names that are in all caps and put them in lowercase
    if ($$name =~ /[A-Z]{4,}/) { $$name =~ tr/A-Z/a-z/ }
    # but uppercase these words
    if ($$name =~ /\b((nadp?h?|atp|v?i+v?|\w))\b/i) { 
	my $s = uc($1);
	$$name =~ s/$1/$s/;
    }
	
    # capitalize the first letter of protein names
    if ($$name =~ /\b([a-z]{3}([A-Z]|\d+))\b/) {
 	my $gene = $1;
	my $new = ucfirst($gene);
	$$name =~ s/$gene/$new/;
    }

    # get rid of bad text
    $$name =~ s/\[imported\]//;
    $$name =~ s/\"//g;
    $$name =~ s/\([^\)]*$//;
    $$name =~ s/\[[^\]]*$//;
    $$name =~ s/\{[^\}]*$//;

    # look for non-homology based info
    $$name =~ s/(authentic )?frameshift//i;
    $$name =~ s/(authentic )?point mutation//i;
    $$name =~ s/ (selenoprotein|selenocysteine\-? ?containing)//i;
    $$name =~ s/, ?(truncation|degenerate|pseudogene)//i;
    $$name =~ s/, ?interruption\-[nc]//i;
    $$name =~ s/(probable|putative|predicted|hypothetical) (lipoprotein|transmembrane protein)(, ?putative)?//i;
    $$name =~ s/similar to .*//i;
    $$name =~ s/alternate start at bp \d+\;//i;
    $$name =~ s/\bpossible\b/putative/;
    $$name =~ s/\bprecursor\b//i;
    $$name =~ s/\bpartial\b//i;
    $$name =~ s/(chloroplast|mitochondrial)//i;
    $$name =~ s/((n|c)-)?terminal//;

    # Handle spacing issues
    $$name =~ s/\s+/ /g;
    $$name =~ s/^\s+//;
    $$name =~ s/\s+$//;
    $$name =~ s/\s*\.$//;

    # look for bad names
    if ($$name =~ /^orf/i ||
	$$name =~ /^\s*protein\s*$/i ||
	$$name =~ /^\s*putative\s*$/i ||
	$$name =~ /hypothetical[0-9\.kda ]* protein/i ||
	$$name =~ /^\s*hypothetical\s*$/i ||
	$$name =~ /^unknown( protein)?/i ||
	$$name =~ /^\s*$/ || 
	$$name =~ /\(\).*/ || 
	$$name =~ /^\s*\d+\s*$/ || 
	$$name =~ /^AGR/ ||
	$$name =~ /^CG(\d+) gene product/ ||
	$$name =~ /^RX/ ||
	$$name =~ /unnamed protein product/i ||
	$$name =~ /No significant matches/i ||
	$$name =~ /conserved domain protein/i ||
	$$name =~ /\bpredicted\b/i ||
	$$name =~ /\bstructure\b/i ||
	$$name =~ /\bCOG\d+(\:)?/ ||
	$$name =~ /incomplete/i ||
	$$name =~ /\-like/i ||
	$$name =~ /uncharacteri(s|z)ed/i ||
	$$name =~ /\bsimilar\b/i ||
	$$name =~ /\bidentical\b/i
	) {
	$$name = "conserved hypothetical protein";
    }

    if ($$name =~ /\bhomolog(ue)?\b/ && $$name !~ /putative/) { $$name = "putative " . $$name }

    # if name is really long for some reason, truncate it
    if (length($$name)>255) {
	my $trunc_name = substr($$name, 0, 255);
	$$name = $trunc_name;
	warn "Had to truncate $ref->{acc} name in clean_annotation:\n$$name\n";
	print LOG "Had to truncate name in clean_annotation\n";
    }
    # Is there an ec hiding in the name?
    if ($$name =~ s/(\d+(\.(\d+|-)){0,2}((\.\d+){3}|\.-))//) {
	if (!$$ec) { $$ec = $1 }
    }
    $$name =~ s/(\[|\()ec\:?\s*(\]|\))//i; # this could be a leftover from EC finding

    # get rid of initial non-word characters
    $$name =~ s/^[^\w\d\[\(\{\<]*//;

    ## remove leading, trailing spaces:
    $$name =~ s/^ *//;
    $$name =~ s/\,*\.* *$//;

    # try to lowercase first letter
    if ($$name =~ /^[A-Z][a-z]{3,}/) {
	$$name = "\l$$name";
    }

    # if the data doesn't look good, delete it
    if ($$name !~ /\w+/ || length($$name) < 3) { $$name = "" }


    print LOG ".....'$$name\n";

    if ($$gene) {
	print LOG "Cleaning '$$gene'.....\n";
	#leading and trailing spaces
	$$gene =~ s/^ *(.*) *$/$1/;
	# TIGR indexing
	$$gene =~ s/\-?\d+$//;
	# screen for bad data
	if ($$gene !~ /\w+/) { $$gene = "" }
	if (length($$gene) < 3) { $$gene = "" }
	print LOG "......'$$gene'\n";
    }

    if ($$ec) {
	print LOG "Cleaning '$$ec'......\n";
	$$ec =~ s/^ *(.*) *$/$1/;
	if ($$ec =~ /^[^\d\.\-]+$/) {
	    $$ec = "";
	}
	print LOG "......'$$ec'\n";
    }
}

sub compare_info {
    my ($ref, $tref) = @_;
    if (! $ref->{product}) {
	print "Initiating annotation: $tref->{product} $tref->{gene} $tref->{EC_number} from $tref->{acc}\n" if ($VERBOSE);
	foreach my $key(keys %$tref) { $ref->{$key} = $tref->{$key} }
    } else {
	if (
	    # hits with rank < 3 are either experimental or highly curated, so info is to be trusted.
	    # But we don't want conserved hypo unless necessary, so don't substitute for that.
	    ($tref->{rank} < 3 && $tref->{rank} < $ref->{rank} ) ||
	    ($ref->{product} =~ /hypothetical|probable|putative|related|unknown/i &&
	     $tref->{product} !~ /hypothetical|probable|putative|related|unknown/i) &&
	    $tref->{product} !~ /hypothetical protein/ ) {
	    print "Upgrading ($ref->{rank}, $tref->{rank}) annotation to $tref->{product} $tref->{gene} $tref->{EC_number} from $tref->{acc} because of rank or name\n" if ($VERBOSE);
	    foreach my $key(keys %$tref) { $ref->{$key} = $tref->{$key} }
	} elsif ($ref->{rank} == $tref->{rank}) {
	    my ($refc, $trefc);
	    if ($ref->{gene}) { $refc++ }
	    if ($ref->{EC_number}) { $refc++ }
	    if ($tref->{gene}) { $trefc++ }
	    if ($tref->{EC_number}) { $trefc++ }
	    if ($trefc > $refc) {
		print "Upgrading annotation to $tref->{product} $tref->{gene} $tref->{EC_number} from $tref->{acc} because of more annotation ($refc, $trefc)\n" if ($VERBOSE);
		foreach my $key(keys %$tref) { $ref->{$key} = $tref->{$key} }
	    }
	}
    }
}

sub add_annotation_level {
    my($product, $per_id, $length) = @_;
    return($product) if ($product =~ /hypothetical|putative|-related|homolog/i);
##
## want to add the words "putative" or "-related" to those products
## that warrant it. 
##
## "high confidence" if the match is greater than $per_id_cutoff spanning
## greater than $length_cutoff
##
## "putative" if the %identity is less than $per_id_cutoff (35) spanning
## greater than $length_cutoff (80).
##
## "-related" if the match is spanning less than $length_cutoff & over
## $per_id_cutoff.
##
## "low confidence" if the match is less than $per_id_cutoff spanning
## under $length_cutoff
##
    print "$per_id, $per_id_cutoff, $length, $length_cutoff\n" if ($DEBUG);
    if (($per_id < $per_id_cutoff) && ($length >= $length_cutoff)) {
	$product .= ", putative";
    } elsif (($per_id >= $per_id_cutoff) && ($length < $length_cutoff)) {
	$product = "conserved domain protein";
    } elsif ($per_id < $per_id_cutoff && $length < $length_cutoff) {
	$product = "hypothetical protein";
    }
    return($product);
}

sub checkForExperimental {
    my($dbproc, $acc) = @_;
    my($ret, $query, $x, @parts);

    (@parts) = split(/\|/,$acc);
    $acc = $parts[0] . ":" . $parts[1];

    $query = "select count(*) from egad.characterized where accession = \"$acc\" and status_vocab_id != 57";
    ($x) = $dbproc->selectrow_array($query);
    if ($x>0) { 
	print "$acc is has experimental row.\n" if ($VERBOSE);
	return(1);
    }
    return 0;
}
