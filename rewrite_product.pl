#!/usr/bin/perl

#Edit WCN and JMo November 4, 2016 to look for translation coordinates for the weird genes

use strict;
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Long;
use DBI;
use Bio::Factory::FTLocationFactory;
use Bio::SeqFeature::Generic;
use Bio::PrimarySeqI;

#my %arg;
#&getopts('D:p:s:S::n:a:', \%arg);
my $db;
my $user = $ENV{'USER'};
my $password;
my @seq_id;
my $seq_acc;
my $seq_file;
my $set_name;
my $feature_file;
GetOptions("db=s" => \$db,
	   "user=s" => \$user,
	   "password=s" => \$password,
	   "seq_id=i" => \@seq_id,
	   "seq_acc=s" => \$seq_acc,
	   "set_name=s" => \$set_name,
	   "feature_file=s" => \$feature_file,
	   "seq_file=s" => \$seq_file
    );
if (!($db && $password)) { die "
USAGE: rewrite_product.pl -db db -password dbpswd [ -seq_id seq_id -seq_acc seq_acc -set_name set_name -feature_file file_of_feature_ids_or_accs -seq_file file_of_seq_ids_or_accs ]\n"; }

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);

if ($seq_acc) { push @seq_id, &get_seq_id_from_acc($dbh, $seq_acc) }
if ($set_name) { 
    my $set_id = set_name_to_set_id($dbh, $set_name);
    push @seq_id, @{&set_id_to_seq_ids($dbh, $set_id)};
}
if ($seq_file) { push @seq_id,  &read_file($seq_file) }

my $SEQ = &get_sequence_by_seq_id($dbh, @seq_id);
my @fids;
foreach my $seq_id(@seq_id) {
#    $SEQ{$seq_id} = &get_sequence_by_seq_id($dbh, $seq_id);
#    my $seqobj = Bio::Seq->new(-seq => $seqstruct->{$seq_id}->{sequence},
#			    -alphabet => 'dna');
#    $SEQ{$seq_id} = $seqobj;
    my $fids = &get_feature_ids_by_seq_id($dbh, $seq_id);
    push @fids, @$fids;
}

if ($feature_file) { push @fids, &read_file($feature_file) }

my $feat = get_features_by_feature_id($dbh, @fids);
my $locfact = Bio::Factory::FTLocationFactory->new();
foreach my $fid(keys %$feat) {
    if ($feat->{$fid}->{'feat_type'} ne "CDS") { next }
    # Cut gene sequence from contig
    my ($set_id) = sort {$b<=>$a} keys %{$feat->{$fid}->{'location'}};
    my $loc = $feat->{$fid}->{'location'}->{$set_id};
    my $seq_id = $loc->{seq_id};
    if (! defined $SEQ->{$seq_id}) {
	my $seqstruct = &get_sequence_by_seq_id($dbh, $seq_id);
	$SEQ->{$seq_id} = $seqstruct->{$seq_id};
    }
    my $seqobj = $SEQ->{$seq_id};

    # in the feature object in ->{'location'}->{$seq_id}->{'translation_coords'} we have the translation coords.
    # How to feed this to bioperl to get correct translation for pseudogenes/weird genes.
    # Decision point for using translation_coords vs feat_min/feat_max
    if (defined $feat->{$fid}->{'location'}->{$seq_id}->{'translation_coords'} &&
	$feat->{$fid}->{'location'}->{$seq_id}->{'translation_coords'} =~ /\.\./) {
	my $location = $locfact->from_string($feat->{$fid}->{'location'}->{$seq_id}->{'translation_coords'});
	my $featobj = Bio::SeqFeature::Generic->direct_new;
	$featobj->location($location);
	$featobj->attach_seq($seqobj);
	$featobj->add_tag_value('condon_start',$feat->{$fid}->{'location'}->{$seq_id}->{'phase'}+1);
	my $prot = $featobj->spliced_seq->translate->seq;
	$prot=~s/\*$//;
    my $q = "UPDATE sequence_features set product=\"$prot\" where feature_id=$fid";
    $dbh->do($q);
	next;
    }
	else { my $subseqobj = $seqobj->trunc($loc->{'feat_min'},
				       $loc->{'feat_max'});
	my $rev;
	if (defined $loc->{'strand'}) {
	    $rev = $loc->{'strand'} < 0 ? 1 : 0;
	}
	my $featseqobj = $rev ? $subseqobj->revcom : $subseqobj;
	
	if ($featseqobj->seq !~ /^[atg]tg/i) {
	    warn "$seq_id $fid does not have usual start: " . $featseqobj->seq . "\n";
	    next;
	}
	if ($featseqobj->seq !~ /(t|u)(aa|ag|ga)$/i) {
	    warn "$seq_id $fid does not have usual stop: " . $featseqobj->seq . "\n";
	    next;
	}
	# Translate sequences to aa if applicable
	my $complete = $loc->{'min_partial'} ||
	    $loc->{'max_partial'} ? 0 : 1; 
	my $protobj = $featseqobj->translate(-complete => $complete,
					     -frame    => $loc->{'phase'});
	my $prot = $protobj->seq;
	if ($prot =~ /\*./) {
	    warn "$seq_id $fid has internal stops: $prot\n";
	    print ">$fid\n" . $feat->{$fid}->{'product'} . "\n" . $prot . "\n";
	    next;
	}
	my $q = "UPDATE sequence_features set product=\"$prot\" where feature_id=$fid";
	$dbh->do($q);
    }
}
exit();

sub read_file {
    my $file = shift;
    open (my $fh, $file) or die "Can't read file '$file': $!\n";

    my @array;
    while (my $line = <$fh>) {
	chomp $line;
	my @f = split(/\s+/, $line);
	push @array, $f[0] if (defined $f[0]);
    }
    return @array;
}
