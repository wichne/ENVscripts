#!/usr/bin/perl
use strict;
use lib $ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;
use Cwd;

my $arg = {};
&getopts('D:i:n:o:', $arg);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg->{D}", 'access', 'access');

my $setid = $arg->{i};
my $setname = $arg->{n};

my $path = ""; # cwd();
my %pt = ("CDS" => "P",
	  "tRNA" => "TRNA",
	  "rRNA" => "RRNA",
	  "tmRNA" => "MISC-RNA",
	  "ncRNA" => "MISC-RNA");

if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $seqref = &get_sequences_by_set_id($dbh, $setid);

open(ELEM, ">genetic-elements.dat");
print ELEM "ID\tCHROM\nNAME\tChromosome\nTYPE\t:CHRSM\nCIRCULAR?\tN\n";
my $elem_txt;

foreach my $seq_id (sort {$a<=>$b} keys %$seqref) {
    my $name;
    foreach my $a(@{$seqref->{$seq_id}->{'accessions'}}) {
	if ($a =~ /PNNL/) { $name = $a }
	else { $name = $a if (! $name) }
    }
    $name =~ s/.*\|//;
    if (! $name) { die "Couldn't get name from accessions: " 
		       . join(" ", @{$seqref->{$seq_id}->{'accessions'}}) . "\n" }

    my $ann_file = $name . ".pf";
#    my $seq_file = $name . ".fsa";
    open (my $annout, ">$ann_file");
#    my $faout = Bio::SeqIO->new(-file => ">$seq_file", -format => 'fasta');
#    $faout->write_seq($mol->{seqobj});

    print ELEM "CONTIG\t" . $name . "\n";
    $elem_txt .= "ID\t" . $name . "\nNAME\t" .$name . "\nTYPE\t:CONTIG\n";
    $elem_txt .= "ANNOT-FILE\t$ann_file\n";
#    $elem_txt .= "SEQ-FILE\t$path/$seq_file\n";
    $elem_txt .= "//\n";

    my $featids = &get_feature_ids_by_seq_id($dbh, $seq_id);
    my $protref = &get_features_by_feature_id($dbh, @$featids);
    foreach my $fid (sort {$a<=>$b} keys %$protref) {
	# make the identifier the locus_tag
	my $acc = $protref->{$fid}->{'accessions'}->{'locus_tag'}->[0];
#	foreach my $src(keys %{$protref->{$fid}->{'accessions'}}) {
#	foreach my $a(@{$protref->{$fid}->{'accessions'}}) {
#	    if ($a =~ /IMG/) { $acc = $a }
#	    else { $acc = $a if (! $acc) }
#	}
	if (! $acc) { warn "No locus_tag for $fid : $acc!!!\n";}

	if (length($protref->{$fid}->{'product'}) == 0 ||
	    !defined($protref->{$fid}->{'product'})) {
	    warn "No product string for feature $fid ($acc). Skipping...";
	    next;
	}
	my ($product, $gene_sym, $ec);
	my ($rank) = sort {$a <=> $b} keys %{$protref->{$fid}->{'annotation'}};
	my ($src) = sort {$a cmp $b} keys %{$protref->{$fid}->{'annotation'}->{$rank}};
	while (my ($k,$v) = each %{$protref->{$fid}->{'annotation'}->{$rank}->{$src}}) {
	    if ($k eq "product") { $product = $v->[0] . " " . $product; }
	    if ($k eq "gene_sym") { $gene_sym = $v->[0] };
	    if ($k eq "EC") { push @$ec, @$v; }
	}

	
	my $loc = $protref->{$fid}->{'location'}->{$setid};
	my ($end5, $end3) = $loc->{'strand'} == 1 
	    ? ($loc->{'feat_min'}, $loc->{'feat_max'}) 
	    : ($loc->{'feat_max'}, $loc->{'feat_min'});
	
	print $annout "ID\t" . $acc . "\n";
	print $annout "NAME\t" . $gene_sym . "\n";
	print $annout "STARTBASE\t" . $end5 . "\n";
	print $annout "ENDBASE\t" . $end3 . "\n";
	print $annout "FUNCTION\t" . $product . "\n";
	print $annout "PRODUCT-TYPE\t$pt{$protref->{$fid}->{feat_type}}\n";
	if (defined $ec) {
	    foreach my $e (@{$ec}) {
		print $annout "EC\t$e\n";
	    }
	}
    print $annout "//\n";
    }
}
print ELEM "//\n" . $elem_txt;
close ELEM;
