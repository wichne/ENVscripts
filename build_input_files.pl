#!/usr/bin/perl

use strict;
use DBI;
use Bio::Seq;
use Getopt::Std;
use Bio::SeqIO;

my $arg = {};
&getopts('D:u:p:s:SIi:', $arg);
my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$arg->{D}", $arg->{u}, $arg->{p}) or die "Check input db credentials (-D -u -p).\n";

print STDERR "Using $arg->{D} as database\n";
print STDERR "Using asmbl_id $arg->{s}\n" if (defined $arg->{s});
if (defined $arg->{S}) {
    print STDERR "Writing nucleotide sequence file (-S)\n";
} else { print STDERR "Not writing nucleotide sequence file (-S)\n"; }
if (defined $arg->{I}) {
    print STDERR "Writing pf info file (-I)\n";
} else { print STDERR "Not writing pf info file (-I)\n"; }
print STDERR "Using index $arg->{i}\n" if (defined $arg->{i});

my $output = $arg->{s} ? $arg->{D} . "_" . $arg->{s} . ".pep" : $arg->{D} . ".pep";
my $seq_output = $arg->{s} ? $arg->{D} . "_" . $arg->{s} . ".seq" : $arg->{D} . ".seq";
my $info_output = $arg->{s} ? $arg->{D} . "_" . $arg->{s} . ".info" : $arg->{D} . ".info";
my $index = "S" . $arg->{i};

my $outr = Bio::SeqIO->new(-file => ">$output",
			   -format => 'fasta',
			   -width => 10000);

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

my ($seqout, $info_file);
$seqout = Bio::SeqIO->new(-file => ">$seq_output",
			  -format => 'fasta',
			  -width => 30000) if ($arg->{S});
open $info_file, ">$info_output";

my @seq_id;
if (! $arg->{s} ) {
    my $seq_id_q = "SELECT asmbl_id FROM assembly ORDER BY asmbl_id";
    my $seq_id_r = $dbh->selectcol_arrayref($seq_id_q);
    @seq_id = @$seq_id_r;
} else { @seq_id = ($arg->{s} ) }


foreach my $seq_id (@seq_id) {
    my $seq_q = "SELECT sequence FROM assembly WHERE asmbl_id = $seq_id";
    my ($seq) = $dbh->selectrow_array($seq_q);
    my $seqo = Bio::Seq->new(-display_id => $seq_id,
			     -seq => $seq);
    
    my $feature_q = "SELECT f.feat_name, end5, end3, locus, product, ec"
	. " FROM feature f, ident i"
	. " WHERE asmbl_id=$seq_id"
	. " AND feat_type='CDS'"
	. " AND f.feat_name=i.feat_name"
	;
    my $feature_ref = $dbh->selectall_hashref($feature_q, 'feat_name');
    my $trans_q = "SELECT f.feat_name, coords, frame, exception"
	. " FROM feature f, translation t"
	. " WHERE asmbl_id=$seq_id"
	. " AND feat_type='CDS'"
	. " AND f.feat_name=t.feat_name";
    my $trans_ref = $dbh->selectall_hashref($feature_q, 'feat_name');
    foreach my $fn (keys %$trans_ref) {
	foreach my $f (keys %{$trans_ref->{$fn}}) {
	    $feature_ref->{$fn}->{$f} = $trans_ref->{$fn}->{$f};
	}
    }
    
    foreach my $ref (sort {$a->{end5} <=>$b->{end5}} values %$feature_ref) {
	if (! $ref->{'locus'}) {
#	    warn "No locus for $ref->{feat_name}. Using ", $arg->{D} , "_$ref->{feat_name}.";
	    $ref->{'locus'} = $arg->{D} . "_" . $ref->{'feat_name'}; }
	if (! defined $ref->{'frame'}) { $ref->{'frame'} = 1 }

	my $seq = $ref->{'end5'} < $ref->{'end3'} ?
	    $seqo->trunc($ref->{'end5'}, $ref->{'end3'}) :
	    $seqo->trunc($ref->{'end3'}, $ref->{'end5'})->revcom;

	my $transl_table = 11;
#	my $description = $ref->{'product'};
	my @ec = split/\s+/, $ref->{'ec'};
#	if (@ec) { $description .= " [EC_number=$_]" foreach (@ec) }

	# This gets the precise na sequence that translates to any 'repaired' pseudogene.
	my $cdsseq = $seq->seq;
	if (defined $ref->{'coords'}) {
	    foreach my $seg(split(/;/,$ref->{'coords'})) {
		my ($end5, $end3) = split/\.{2}/, $seg;
		if ($end5 < $end3) {
		    my $segseq = $seqo->trunc($end5, $end3);
		    $cdsseq .= $segseq->seq;
		} else {
		    my $segseq = $seqo->trunc($end3, $end5)->revcom;
		    $cdsseq = $segseq->seq . $cdsseq;
		}
	    }
	}

	my $locus = defined $arg->{'i'} ? "${index}_" . $ref->{'locus'} : $ref->{'locus'};
	my $geneObj = Bio::Seq->new(-display_id => $locus,
#				    -description => $description,
				    -seq => $cdsseq);

	# Now we need to check the geneObj sequence for internal stops. I know, they're not supposed to be there, but they are and it screws up MSOAR downstream becuase mafft screens out the termination characters which makes an inconsistency between the aa and nuc sequence in pal2nal.
	my $checkseq = $geneObj->seq;
	my $finalseq;
	while ($checkseq =~ /(.{3})/g) {
	    my $codon = $1;
	    if ($codon !~ /t(a(a|g))|(ga)/i) { $finalseq .= $codon }
	}
	$geneObj->seq($finalseq);

	my $protObj = $geneObj->translate(-frame => $ref->{'frame'} - 1,
					  #-complete => 1,
					  -codontable_id => $transl_table); # codontable_id

	if ($ref->{'exception'}) {
	    my @te = split/\;/, $ref->{'exception'};
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
	$outr->write_seq($protObj);

	if ($arg->{S}) {
#	    $seq->display_id($locus);
	    $seqout->write_seq($geneObj);
	}

	if ($arg->{I}) {
	    my $dir = $ref->{'end5'} < $ref->{'end3'} ? "+" : "-";
	    printf $info_file "%s\t%s\t%s\t%s\t%s\n",
	    ($locus, $locus, $seq_id, $dir, $ref->{'end5'});
	}
    }
}

exit();
