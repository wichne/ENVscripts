#!/usr/bin/perl
#JMo 09/24/2015, change lib to /scripts, updated subroutine (line 85) to reflect ENV.pm

use lib "/scripts";
use ENV;
use Getopt::Std;
use strict;

my $DEBUG = 0;

my %arg;
&getopts("D:p:u:g:f:i:n:h", \%arg);

my $dbh = &connect(\%arg);

my $genbank_file = $arg{g};
my $seq_file;
my $prefix;
if ($genbank_file =~ /(.*)\.gb[fk]?/) {
    $prefix = $1;
    $seq_file = $prefix . ".seq";
} else {
    die "$genbank_file has wrong suffix for a genbank file (.gb, .gbf, or .gbk)\n";
}

# make the .seq file
&run_command_to_make_file($seq_file, "gb2fasta.pl $genbank_file");

# Get the sequences of the existing features
my $old_fa;
my $set_name;
if ($arg{'i'}) {
    $set_name = &set_id_to_name($dbh, $arg{i});
    $old_fa = $set_name . ".seq";
    &run_command_to_make_file($old_fa, "/share/scripts/devel/ENV/set_to_seq.pl -D $arg{D} -i $arg{i} -o $old_fa");
} elsif ($arg{'n'}) {
    $set_name = $arg{'n'};
    $old_fa = $arg{'n'} . ".seq";
    &run_command_to_make_file($old_fa, "/share/scripts/devel/ENV/set_to_seq.pl -D $arg{D} -n $arg{n} -o $old_fa");
} elsif ($arg{'f'}) {
    $old_fa = $arg{'f'};
    if (!(-e $old_fa && -r $old_fa)) { die "Can't open fastafile $old_fa: $!\n"; }
}

# Run nucmer
my $delta_file = $prefix . "_v_" . $set_name;
&run_command_to_make_file("$delta_file.delta", "nucmer -maxmatch -p $delta_file $seq_file $old_fa");

# use the show-coords output to do an initial mapping
open my $NUCMER, "show-coords -lrHT $delta_file.delta | ";
my $perfect;
my %MAP;
my %REV;
while (my $l=<$NUCMER>) {
    chomp $l;
    my ($seq_lo,
	$seq_hi,
	$feat_start,
	$feat_end,
	$seq_reg_len,
	$feat_reg_len,
	$perid,
	$seq_len,
	$feat_len,
	$ref_acc,
	$feat_acc) = split/\t/, $l;
    $MAP{$ref_acc}->{$feat_acc} = { 'perid' => $perid/100,
				    'perlen' => $feat_reg_len/$feat_len };
    $REV{$feat_acc}->{$ref_acc} = { 'perid' => $perid/100,
				    'perlen' => $feat_reg_len/$feat_len };
    if ($feat_reg_len/$feat_len == 1 &&
	$perid == 100) {
	$perfect++;
    }
}

my $gbo = Bio::SeqIO->new(-file => $genbank_file,
			  -format => 'genbank');
while (my $seqo = $gbo->next_seq) {
    my @features = $seqo->get_SeqFeatures();
    # the first feature in a genbank file is the sequence itself
    # we need to insert the sequence
    my $seqObj = shift @features;
    # check to see if the sequence is already loaded
    my $seq_id = &get_seq_id_by_seq_accession($dbh, $seqObj->seq->display_name);
    if (! $seq_id) {
	$seq_id = &load_sequence_SeqFeature($dbh, $seqObj);
    }
    print STDERR "Sequence " . $seqObj->display_name . " is seq_id: $seq_id\n";

    my ($start, $end, $locus_tag);
    foreach my $featObj (@features) {
	if ($featObj->primary_tag eq "sig_peptide" ||
	    $featObj->primary_tag eq "repeat_region") { next }
	if ($featObj->primary_tag eq "gene") {
	    if ($featObj->has_tag('locus_tag')) {
		$locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
		$start = $featObj->start;
		$end = $featObj->end;
	    }
	    next;
	}
	
	if (!($featObj->has_tag("locus_tag"))) {
	    if ($featObj->start eq $start && $featObj->end eq $end) {
		$featObj->set_attributes(-tag => {'locus_tag' => $locus_tag });
	    }
	} else {
	    $locus_tag = [$featObj->get_tag_values('locus_tag')]->[0];
	}

	my $fid = &get_feature_id_by_accession($dbh, $locus_tag);
	if ($fid) {
	    print STDERR "$locus_tag maps to $fid\n";
	    &load_seq_feat_mappings($dbh, $fid, $seq_id, $featObj) unless ($DEBUG);
	    $locus_tag = "";
	} elsif (defined $MAP{$locus_tag}) {
	    foreach my $facc(sort { $MAP{$locus_tag}->{$b}->{'perid'} * $MAP{$locus_tag}->{$b}->{'perlen'} <=>
					$MAP{$locus_tag}->{$a}->{'perid'} * $MAP{$locus_tag}->{$a}->{'perlen'} }
			     keys %{$MAP{$locus_tag}}) { 

		# how good is this hit?
		my $ref = $MAP{$locus_tag}->{$facc};
		if (!($ref->{'perid'} == 1 &&
		    $ref->{'perlen'} == 1 &&
		    scalar(keys %{$REV{$facc}}) == 1)) {
		    printf "Best hit to %s %s is %s (which has %d other hits) with %.2f identity across %.2f of length. Should this mapping occur? (Y/n) ", ($featObj->primary_tag, $locus_tag, $facc, scalar(keys %{$REV{$facc}})-1,  $ref->{perid} ,$ref->{perlen});
		    my $answer = <STDIN>;
		    if ($answer =~ /n/i) { next }
		}
		my @acc = split/\|/, $facc;
		foreach my $x(@acc) { 
		    next if ($x =~ /gnl|fig/);
		    $fid = &get_feature_id_by_accession($dbh, $x);
		    if ($fid) { last }
		}
		if (! $fid) { 
		    print STDERR "Why can't I get a feature_id from $facc?";
		} else { 
		    print STDERR "$locus_tag maps to $fid\n";
		    last;
		}
	    }

	    if ($fid) {
		&load_seq_feat_mappings($dbh, $fid, $seq_id, $featObj) unless ($DEBUG);
	    } else {
		print STDERR "Loading $locus_tag as new gene\n";
		my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj) unless ($DEBUG);
		&load_feature_annotations($dbh, $feat_id, $featObj->annotation) unless ($DEBUG);
	    }
	    $locus_tag = "";
	} else {
	    if (! $locus_tag) {
		printf "Found %s at %d/%d. Should I load it as a new feature? (y/N) ", ($featObj->primary_tag, $featObj->start, $featObj->end);
		my $ans = <STDIN>;
		chomp $ans;
		if (! $ans || $ans =~ /n/i) { next }
	    }
	    print STDERR "Loading $locus_tag as new gene\n";
	    my $feat_id = &load_SeqFeature($dbh, $seq_id, $featObj, $seqObj) unless ($DEBUG);
	    &load_feature_annotations($dbh, $feat_id, $featObj->annotation) unless ($DEBUG);
	    $locus_tag = "";
	}
    }
}


sub run_command_to_make_file {
    my $file = shift;
    my $cmd = shift;
    if (-e $file) {
	print "File $file already exists. Remake it (y/N) ";
	my $ans = <STDIN>;
	chomp $ans;
	if ($ans =~ /y/i) {
	    system($cmd);
	}
    } else {
	system($cmd);
    }
}
