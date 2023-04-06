#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use lib $ENV{ENVSCRIPTS};
use SGDB;
use strict;
use Getopt::Std;

my $args = {};
&getopts('D:u:p:i:F:a:', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $asmbl_id = $args->{'a'};

my %AA = ('Ala' => 'A',
	  'Cys' => 'C',
	  'Asp' => 'D',
	  'Glu' => 'E',
	  'Phe' => 'F',
	  'Gly' => 'G',
	  'His' => 'H',
	  'Iso' => 'I',
	  'Lys' => 'K',
	  'Leu' => 'L',
	  'Met' => 'M',
	  'Asn' => 'N',
	  'Pyl' => 'O',
	  'Pro' => 'P',
	  'Gln' => 'Q',
	  'Arg' => 'R',
	  'Ser' => 'S',
	  'Thr' => 'T',
	  'Sel' => 'U',
	  'Val' => 'V',
	  'Trp' => 'W',
	  'Tyr' => 'Y');

my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);

my $dbh = DBI->connect("dbi:mysql:server=localhost;database=$db",
		       $user, $pswd, {'AutoCommit'=> 0});

while (my $seqo = $in->next_seq) {
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence into assembly, stan and asmbl_data
    my $seqf = shift @features;
    if (! $asmbl_id) { $asmbl_id = &load_assembly($seqf); }
    foreach my $featObj (@features) {
	my $feat_type = $featObj->primary_tag;
	if ($feat_type eq "CDS") {
	    my @coords;
	    if ( $featObj->location->isa('Bio::Location::SplitLocationI') ) {
		for my $location ( $featObj->location->sub_Location ) {
		    if ($location->strand > 0) {
			push @coords, join("..",($location->start, $location->end));
		    } else {
			push @coords, join("..",($location->end, $location->start));
		    }
		}
	    } else {
		@coords = $featObj->strand > 0 ?
		    (join("..", ($featObj->start, $featObj->end))) :
		    (join("..", ($featObj->end, $featObj->start)));
	    }

	    my %transl = ('coords' => join(";",@coords),
			  'frame' => 1);

	    my $target_tx;
	    my @tags = $featObj->get_all_tags;
# codon_start - [1,2,3] single value
# exception - for RNA editing notes; multiple values?
# pseudogene - ['processed', 'unprocessed', 'unitary', 'allelic', 'unknown']  one value
	    # processed - reverse tx of mRNA and re-integration into genome
	    # unprocessed - degraded copy of duplicated gene
	    # unitary - degraded copy of gene
	    # allelic - functional alleles exist in population
	    # unknown - 
# ribosomal_slippage - no value
# translation - the translation (use it to check info)	one value
# transl_except - pos:<location>,aa:<amino_acid>; multiple values
# transl_table - number; one value
# trans_splicing - no value
	    foreach my $tag (@tags) {
		my @values = $featObj->get_tag_values($tag);
		if ($tag eq "codon_start") { $transl{'frame'} = $values[0]; }
		elsif ($tag eq "locus_tag") { $transl{'feat_name'} = &get_feat_name_from_locus($dbh, $values[0]) }
		elsif ($tag eq "exception") {
		    foreach my $v(@values) {
			$transl{'comment'} .= "$v;";
		    }
		} elsif ($tag =~ /^pseudo/) { $transl{'qualifier'} .= "$tag=$values[0];"}
		elsif ($tag eq "ribosomal_slippage") { $transl{'qualifier'} .= "$tag;" }
		elsif ($tag eq "translation") { $target_tx = $values[0] }
		elsif ($tag eq "transl_except") {
		    foreach my $v(@values) {
			$transl{'exception'} .= "$v;" }
		} elsif ($tag eq "transl_splicing") { $transl{'qualifier'} .= "$tag;" }
	    }
	    if (! defined $transl{'feat_name'}) { next; }
	    # Now try out the info to see if we get back the reported translation
	    my $featseq;
	    my $protseq;
	    foreach my $seg(split(/;/,$transl{'coords'})) {
		my ($end5, $end3) = split/\.{2}/, $seg;
		if ($end5 < $end3) {
		    my $segseq = $seqo->trunc($end5, $end3);
		    $featseq .= $segseq->seq;
		} else {
		    my $segseq = $seqo->trunc($end3, $end5)->revcom;
		    $featseq = $segseq->seq . $featseq;
		}
	    }
	    my $geneObj = Bio::Seq->new(-display_id => $transl{'feat_name'},
					-seq => $featseq);
	    my $protObj = $geneObj->translate(-frame => $transl{'frame'} - 1); # codontable_id

	    if ($transl{'exception'}) {
		my @te = split/\;/, $transl{'exception'};
		foreach my $te (@te) {
		    if ($te =~ /\:(\d+)(\.{2}(\d+))?\,aa\:(\w+)/) {
			my $lo = $1;
			my $hi = $3;
			my $aa = $4;
			if ($aa eq "TERM") {
			    my $newprot = $protObj->trunc(1,($lo+2)/3 - 1);
			    $protObj->seq($newprot->seq);
			} else {
			    my $pos = ($lo + 2) / 3;
			    my @prot = split(/\s*/,$protObj->seq);
			    $prot[$pos] = $AA{$aa};
			    $protObj->seq(join("",@prot));
			}
		    } else { warn "Unexpected format for transl_except: $te\n"; }
		}
	    }

	    my $prot = $protObj->seq;
	    $prot =~ s/^./M/ if ($geneObj->seq =~ /^[ATG]TG/);
	    $prot =~ s/\*$//;
	    $protObj->seq($prot);

	    if ($transl{'qualifier'} !~ /pseudo/ &&
		($protObj->seq =~ /\*\w/ ||
		(defined $target_tx && $protObj->seq ne $target_tx))) {
		warn "$transl{feat_name} : asmbl_id=$asmbl_id coords=$transl{coords} $transl{exception}\n"
		    . "yields translation:\n" . $protObj->seq . "\nTarget:\n" . $target_tx . "\n\n";
	    }
	}
    }
}
