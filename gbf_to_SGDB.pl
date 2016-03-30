#!/usr/bin/perl

use DBI;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Std;
use Data::Dumper;
use strict;
use lib $ENV{SCRIPTS};
#use lib $ENV{TIGR_SCRIPTS};
use SGDB;
use vars qw ($dbh %FEAT_NAME);

my $args = {};
&getopts('D:u:p:i:F:a:C', $args);

my $filename = $args->{'i'} or die "Need to provide filename with -i\n";
my $format = $args->{'F'} ? $args->{'F'} : undef;
my $db = $args->{'D'} or die "Need to provide database name with -D\n";
my $pswd = $args->{'p'} or die "Need to provide database password with -p\n";
my $user = $args->{'u'} ? $args->{'u'} : $ENV{'user'};
my $asmbl_id = $args->{'a'};

my $in = Bio::SeqIO->new(-file => $filename,
			 -format => $format);

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'bladerunner';
my $user = $ENV{USER};
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", $user, $pswd, {'AutoCommit'=>0});

if ($args->{'C'}) {
    my @tables = ("assembly",
		  "asmbl_data",
		  "stan",
		  "ident",
		  "feature",
		  "translation");

    foreach my $t (@tables) {
	my $q = "TRUNCATE $t";
	$dbh->do($q);
    }
}

my $q = "select locus from ident";
my $existing = $dbh->selectall_hashref($q, 'locus');

while (my $seqo = $in->next_seq) {
    my %FEAT;
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence into assembly, stan and asmbl_data
    my $seqf = shift @features;
    $asmbl_id = &load_assembly($seqf) unless ($args->{'a'});
    print STDERR "Loading " . $seqo->display_id . " -> $asmbl_id...\n";
    foreach my $featObj (@features) {
	&process_feature($asmbl_id, $featObj, \%FEAT);
    }
    foreach my $locus(keys %FEAT) {
	next if (defined $existing->{$locus});
	my $feat_type = $FEAT{$locus}->{'feature'}->{'feat_type'};
	my $feat_name = &get_next_feat_name($dbh, {'feat_type'=>$feat_type});
	$FEAT{$locus}->{'feature'}->{'feat_name'} = $feat_name;
	&insert_feature($dbh, $FEAT{$locus}->{'feature'});
	$FEAT{$locus}->{'ident'}->{'feat_name'} = $feat_name;
	&insert_ident($dbh, $FEAT{$locus}->{'ident'});
	$FEAT{$locus}->{'translation'}->{'feat_name'} = $feat_name;
	&insert_translation($dbh, $FEAT{$locus}->{'translation'});
    }
}

print "\n";
exit();

#-----------> Subroutines

sub load_assembly {
    my $assemblyObj = shift;

    my %assembly = ('sequence' => $assemblyObj->seq->seq);
    my $asmbl_id = &insert_assembly($dbh, \%assembly);

    my %asmbl_data = ('name' => $assemblyObj->seq->description,
		      'type' => $assemblyObj->get_tag_values("mol_type"),
		      'topology' => "linear",
		      'firewall_toggle' => 0,
		      'accession_toggle' => 0);
    my $asmbl_data_id = &insert_asmbl_data($dbh, \%asmbl_data);

    my %stan = ('asmbl_id' => $asmbl_id,
		'asmbl_data_id' => $asmbl_data_id,
		'original_asmbl_id' => $asmbl_id,
		'iscurrent' => 1,
		'change_log' => 0);
    &insert_stan($dbh, \%stan);
    return $asmbl_id;
}

sub process_feature {
    # later we will insert feature, translation and ident, so the purpose of
    # this subroutine is to build those data structures. One critical aspect
    # is combining the gene and CDS features
    my $asmbl_id = shift;
    my $featObj = shift;
    my $FEAT = shift;

    my $feat_type = $featObj->primary_tag;
    my @tags = $featObj->get_all_tags;
    # as of 2013-04-15, valid 'CDS' tags are:
#                       /allele="text"
#                       /artificial_location="[artificial_location_value]"
#                       /citation=[number]
#   feature             /codon_start=<1 or 2 or 3>
#                       /db_xref="<database>:<identifier>"
#   ident               /EC_number="text"
#   translation/prop    /exception="[exception_value]"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#   ident               /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#   ident               /locus_tag="text" (single token)
#                       /map="text"
#   ident               /note="text"
#                       /number=unquoted text (single token)
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#   ident               /product="text"
#   feature_accessions  /protein_id="<identifier>"
#   ident/transl/prop   /pseudo
#   ident/transl/prop   /pseudogene="TYPE"
#   transl/prop         /ribosomal_slippage
#                       /standard_name="text"
#   feature/transl      /translation="text"
#   transl/prop         /transl_except=(pos:<base_range>,aa:<amino_acid>)
#   translation         /transl_table =<integer>
#   transl/prop         /trans_splicing

    # and valid 'gene' tags are:
#                       /allele="text"
#                       /citation=[number]
#                       /db_xref="<database>:<identifier>"
#                       /experiment="[CATEGORY:]text"
#                       /function="text"
#   ident               /gene="text"
#                       /gene_synonym="text"
#                       /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#   ident               /locus_tag="text" (single token)
#                       /map="text"
#   ident               /note="text"
#                       /old_locus_tag="text" (single token)
#                       /operon="text"
#   ident               /product="text"
#   ident/prop          /pseudo
#   ident/prop          /pseudogene="TYPE"
#                       /phenotype="text"
#                       /standard_name="text"
#   transl/prop         /trans_splicing

    # so %FEAT will look like $FEAT{locus_tag}->{'ident'} = { 'field' => 'value' }
    #                                         ->{'feature'} = { 'field' => 'value' }
    #                                         ->{'translation'} = {'field' => 'value'}
    # which means if we don't have a locus_tag, we're in trouble...

    if (! grep /locus_tag/, @tags) {
#!!!	&load_feature ($asmbl_id, $featObj);
	return;
    } else {
	if ($feat_type eq "misc_feature") { next }
	my ($locus) = $featObj->get_tag_values('locus_tag');
	my ($end5, $end3) = $featObj->strand > 0 ?
	    ($featObj->start, $featObj->end) : 
	    ($featObj->end, $featObj->start);

	$feat_type = "CDS" if ($feat_type eq "gene");
	$FEAT->{$locus}->{'feature'}->{'asmbl_id'} = $asmbl_id;
	$FEAT->{$locus}->{'feature'}->{'feat_type'} = $feat_type;
	$FEAT->{$locus}->{'feature'}->{'end5'} = $end5;
	$FEAT->{$locus}->{'feature'}->{'end3'} = $end3;
	$FEAT->{$locus}->{'feature'}->{'sequence'} = $featObj->seq->seq;
	$FEAT->{$locus}->{'feature'}->{'feat_method'} = $featObj->source_tag;
	
	foreach my $tag (@tags) {
	    my @values = $featObj->get_tag_values($tag);
	    if ($tag eq "EC_number") {$FEAT->{$locus}->{'ident'}->{'ec'} = join " ", @values;}
	    elsif ($tag eq "gene") {$FEAT->{$locus}->{'ident'}->{'gene_sym'} = join "; ", @values;}
	    elsif ($tag eq "locus_tag") {$FEAT->{$locus}->{'ident'}->{'locus'} = $values[0];}
	    elsif ($tag eq "note") {$FEAT->{$locus}->{'ident'}->{'comment'} = join "; ", @values;}
	    elsif ($tag eq "product") {
#	    elsif ($value =~ /\b(\d{1,3}(\.[\d,\-]{1,4}){1,3})/) {
		#	$ec .= $1 . " ";
#	    }
#	    elsif ($value =~ /\b(GO:\d{7})\b/) {
#		push @GO, $1;
#	    }
		$FEAT->{$locus}->{'ident'}->{'product'} = $values[0];
	    }
	    elsif ($tag eq "translation") {
		$FEAT->{$locus}->{'feature'}->{'protein'} = $values[0];
		$FEAT->{$locus}->{'translation'}->{'protein'} = $values[0];}
	    elsif ($tag eq "codon_start") { $FEAT->{$locus}->{'translation'}->{'frame'} = $values[0]; }
	    elsif ($tag eq "exception") {
		foreach my $v(@values) {
		    $FEAT->{$locus}->{'translation'}->{'comment'} .= "$v;";
		}
	    }
	    elsif ($tag =~ /^pseudo/) { $FEAT->{$locus}->{'translation'}->{'qualifier'} .= "$tag=$values[0];"}
	    elsif ($tag eq "ribosomal_slippage") { $FEAT->{$locus}->{'translation'}->{'qualifier'} .= "$tag;" }
	    elsif ($tag eq "transl_except") {$FEAT->{$locus}->{'translation'}->{'exception'} = join ";", @values;}
	    elsif ($tag eq "transl_splicing") { $FEAT->{$locus}->{'translation'}->{'qualifier'} .= "$tag;" }
	}

	$FEAT->{$locus}->{'translation'}->{'coords'} = $end5 . ".." . $end3 if (! defined $FEAT->{$locus}->{'translation'}->{'coords'});	    if (! $FEAT->{$locus}->{'translation'}->{'frame'}) { $FEAT->{$locus}->{'translation'}->{'frame'} = 1 }

	if ($featObj->primary_tag eq "CDS") {
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
	    $FEAT->{$locus}->{'translation'}->{'coords'} = join ";", @coords;
	}

	if (! $FEAT->{$locus}->{'ident'}->{'product'}) {
#	    print STDERR "No 'product' for $feat_type '$locus'.";
	    if ($FEAT->{$locus}->{'ident'}->{'comment'}) {
#		print STDERR "\tUsing 'comment' -> $FEAT->{$locus}->{ident}->{comment}\n";
		$FEAT->{$locus}->{'ident'}->{'product'} = $FEAT->{$locus}->{'ident'}->{'comment'};
	    } else {
#		print STDERR "\tUsing 'unknown'.\n";
		$FEAT->{$locus}->{'ident'}->{'product'} = "unknown";
	    }
	}
    }
}

sub load_feature {
    my $asmbl_id = shift;
    my $featObj = shift;
    
    my ($ident, $transl) = &parse_tags($featObj);

    my $feat_type = $featObj->primary_tag;

    my $feat_name = &get_next_feat_name($feat_type);
    

    my ($end5, $end3) = $featObj->strand > 0 ?
	($featObj->start, $featObj->end) : 
	($featObj->end, $featObj->start);

    if (! $ident->{'locus'}) {
#	warn "No locus_tag for $feat_type $end5..$end3\n";
    }

    my %feature = ('asmbl_id' => $asmbl_id,
		   'feat_type' => $feat_type,
		   'feat_name' => $feat_name,
		   'end5' => $end5,
		   'end3' => $end3,
		   'sequence' => $featObj->seq->seq,
		   'feat_method' => $featObj->source_tag,
		   );
    $feature{'protein'} = $transl->{'protein'} if (defined $transl->{'protein'});

    &insert_feature($dbh, \%feature);
    
    if ($feat_type !~ /(CDS)|(RNA)|(gene)/) { return }
    
    if ($feat_type eq "CDS") {
	$transl->{'feat_name'} = $feat_name;
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
	$transl->{'coords'} = join ";", @coords;
	if (! $transl->{'frame'}) { $transl->{'frame'} = 1 }
	&insert_translation($dbh, $transl);
    }

    if (! $ident->{'product'}) {
#	warn "No 'product' for $feat_type '$ident->{locus}'.";
	if ($ident->{'comment'}) {
#	    warn "\tUsing 'comment' -> $ident->{comment}\n";
	    $ident->{'product'} = $ident->{'comment'};
	} else {
#	    warn "\tUsing 'unknown'.\n";
	    $ident->{'product'} = "unknown";
	}
    }
    $ident->{'feat_name'} = $feat_name;
    &insert_ident($dbh, $ident);
#    $sgdb->insert_GO(\%GO);
}

sub parse_tags {
    my $featObj = shift;
    my @tags = $featObj->get_all_tags;

    my %ident;
    my %transl;

    my ($locus, $protein, $product, $gene_sym, $ec, @GO, $comment);

    foreach my $tag (@tags) {
	my @values = $featObj->get_tag_values($tag);

	if ($tag eq "EC_number") {
	    $ident{'ec'} .= join " ", @values;
	}

	if ($tag eq "gene") {
	    $ident{'gene_sym'} = join "; ", @values;
	}
	
	if ($tag eq "locus_tag") {
	    $ident{'locus'} = $values[0];
	}

	if ($tag eq "note") {
	    $ident{'comment'} = join "; ", @values;
	}

	if ($tag eq "product") {
#	    if ($value =~ /\b(\d{1,3}(\.[\d,\-]{1,4}){1,3})/) {
	#	$ec .= $1 . " ";
#	    }
#	    if ($value =~ /\b(GO:\d{7})\b/) {
#		push @GO, $1;
#	    }
	    $ident{'product'} = $values[0];
	}

	if ($tag eq "translation") {
	    $transl{'protein'} = $values[0];
	}

	if ($tag eq "codon_start") { $transl{'frame'} = $values[0]; }

	if ($tag eq "exception") {
	    foreach my $v(@values) {
		$transl{'comment'} .= "$v;";
	    }
	}

	if ($tag =~ /^pseudo/) { $transl{'qualifier'} .= "$tag=$values[0];"}

	if ($tag eq "ribosomal_slippage") { $transl{'qualifier'} .= "$tag;" }

	if ($tag eq "transl_except") {
	    $transl{'exception'} = join ";", @values;
	}

	if ($tag eq "transl_splicing") { $transl{'qualifier'} .= "$tag;" }
    }
    chop $ec if (defined $ec);
    return (\%ident, \%transl);
}

# sub get_next_feat_name {
#     my ($feat_type) = shift;
#     my $feat_name;
#     if (defined $FEAT_NAME{$feat_type}) {
# 	$feat_name = sprintf "$feat_type%05i", (++$FEAT_NAME{$feat_type});
#     } else {
# 	my $query = "SELECT MAX(feat_name) FROM feature"
# 	    . " WHERE feat_type=\"$feat_type\"";
# 	my ($max_feat_name) = $dbh->selectrow_array($query);
# 	if (defined $max_feat_name && $max_feat_name =~ /(\d+)$/) {
# 	    my $index = $1 + 1;
# 	    $feat_name = sprintf "$feat_type%05i", ($index);
# 	    $FEAT_NAME{$feat_type} = $index;
# 	} else {
# 	    $feat_name = $feat_type . "00001";
# 	    $FEAT_NAME{$feat_type} = 1;
# 	}
#     }
#     return $feat_name;
# }
