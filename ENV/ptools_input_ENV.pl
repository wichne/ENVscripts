#!/usr/bin/perl

use DBI;
use strict;
#no strict 'refs';
use Getopt::Long;
use lib $ENV{SCRIPTS};
use ENV;
use Cwd;

my ($db, @id, $wgs, $min_contig_length, $help);
my %param = ("a=i@" => \@id,
	     "d=s"  => \$db,
	     "n=s" => \$set_name,
	     "h" => \$help);

&GetOptions(%param);
my $user     = "access";
my $password = "access";
my $db       = $db;
my $server = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbtype   = "mysql";

my $path = cwd();
open(ELEM, ">genetic-elements.dat");

my $dbh = DBI->connect("dbi:$dbtype:host=$server;db=$db", $user, $password);

#open(ORG, ">organism-params.dat");
my $pid;
until ($pid) {
    print "Enter an ID (2-10 alphanumeric; must start with a letter) for this project (" . $set_name . "): ";
    $pid = <STDIN>;
    chomp $pid;
    $pid = $pid =~ /[a-zA-Z][\w]{1,9}/ ? $pid : undef;
}
my $dbname = lc($pid);
$dbname = "\u$dbname" . "Cyc";

#print ORG "ID\t$pid\n"
#    . "STORAGE\tMySQL\n"
#    . "NAME\t" . $sgc->organism . "\n"
#    . "PRIVATE\tT\n"
#    . "RANK\tspecies\n"
#    . "DOMAIN\tTAX-" . $sgc->taxon_id . "\n"
#    . "CODON-TABLE\t" . $sgc->gcode . "\n"
#    . "DBNAME\t$dbname\n"
#    . "NCBI-TAXON-ID\t" . $sgc->taxon_id . "\n";
#close ORG;

#if (! @id) { 
#    @id = $sgc->current_asmbl_ids
#}

#$sgc->populate_molecules;

my %pt = ("CDS" => "P",
	  "tRNA" => "TRNA",
	  "rRNA" => "RRNA",
	  "tmRNA" => "MISC-RNA",
	  "ncRNA" => "MISC-RNA");

print ELEM "ID\tCHROM\nNAME\tChromosome\nTYPE\t:CHRSM\nCIRCULAR?\tN\n";
my $elem_txt;
foreach my $id (@id) {
    my $mol = $sgc->molecule($id);
    print ELEM "CONTIG\t" . $mol->{name} . "\n";

    my $ann_file = $mol->{name} . ".pf";
    my $seq_file = $mol->{name} . ".fsa";

    $elem_txt .= "ID\t" . $mol->{name} . "\nNAME\t" .$mol->{name} . "\nTYPE\t:CONTIG\n";
    $elem_txt .= "ANNOT-FILE\t$path/$ann_file\n";
    $elem_txt .= "SEQ-FILE\t$path/$seq_file\n";
    $elem_txt .= "//\n";

    my $faout = Bio::SeqIO->new(-file => ">$seq_file", -format => 'fasta');
    $faout->write_seq($mol->{seqobj});

    my $featref = $sgc->populate_all_orfs_on_assembly($id);
    open (my $annout, ">$ann_file");
    foreach my $fn (sort keys %$featref) {
	my $fobj = $featref->{$fn};
	my ($ft, $class) = $fobj->get_feat_type_and_class;
	if (! defined $pt{$ft}) { next }
	my $ann = $fobj->get_annotation;
	print $annout "ID\t" . $fobj->get_locus . "\n";
	print $annout "NAME\t" . $ann->{'gene_sym'} . "\n";
	print $annout "STARTBASE\t" . $fobj->end5 . "\n";
	print $annout "ENDBASE\t" . $fobj->end3 . "\n";
	print $annout "FUNCTION\t" . $ann->{'product'} . "\n";
	print $annout "PRODUCT-TYPE\t$pt{$ft}\n";
	if (defined $ann->{ec}) {
	    foreach my $e (@{$ann->{ec}}) {
		print $annout "EC\t$e\n";
	    }
	}
	print $annout "//\n";
    }
}
print ELEM "//\n" . $elem_txt;
close ELEM;
