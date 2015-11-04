#!/usr/bin/perl
#edit JMo 10022015, add input options

use strict;
use lib "/share/scripts"; #$ENV{SCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use CGI qw/:standard/;
use Bio::SeqIO;

#my $db = &param('db');
#my $setid = &param('set_id');
#my $setname = &param('seq_set_name');
my $opt = {};
&getopts('D:i:n:o:', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = DBI->connect("dbi:mysql:host=$host;db=$db", 'access', 'access');

if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $seq_ids = &set_id_to_seq_ids($dbh, $setid);

my $output;
foreach my $seqid (sort {$a<=>$b} @$seq_ids) {
    my $contig = "\"" . join("\n",@{&get_accession_by_seq_id($dbh, $seqid)}) . "\"";
    my $feat_id_aref = &get_feature_ids_by_seq_id($dbh, $seqid);
    my $feat_ref = &get_features_by_feature_id($dbh, @$feat_id_aref);
    foreach my $featid ( sort { $feat_ref->{$a}->{'location'}->{$seqid}->{'feat_min'} <=>
				    $feat_ref->{$b}->{'location'}->{$seqid}->{'feat_min'} } keys %$feat_ref) {

	my $feat_type = $feat_ref->{$featid}->{'feat_type'};
	# location and size stuff
	my $locref = $feat_ref->{$featid}->{'location'}->{$seqid};
	my ($start, $stop, $strand) = ($locref->{'feat_min'}, $locref->{'feat_max'}, $locref->{'strand'});
	$start = "<" . $start if ($locref->{'min_partial'});
	$stop = $stop . ">" if ($locref->{'max_partial'});
	$strand = $strand < 0 ? "-" : $strand > 0 ? "+" : "";
	my $aa;
	# the length in aa will be either the actual length of the product field, or,
	# if the product field is empty
	if ($feat_ref->{$featid}->{'feat_type'} eq "CDS") {
	    if ($feat_ref->{$featid}->{'product'} =~ /\w+/) {
		$aa = length($feat_ref->{$featid}->{'product'});
	    } else {
		$aa = $locref->{'min_partial'} || $locref->{'max_partial'} ?
		    "(" . sprintf("%.1f",(abs($locref->{'feat_max'} - $locref->{'feat_min'}) + 1)) . "+)" :
		    "(" . sprintf("%.1f",(abs($locref->{'feat_max'} - $locref->{'feat_min'}) + 1)) . ")";
	    }
	}

	# Accessions
	my $PNNL_id = "\"" . join("\n", @{$feat_ref->{$featid}->{'accessions'}->{'PNNL'}}) . "\""
	    if ($feat_ref->{$featid}->{'accessions'}->{'PNNL'});
	my $RAST_id = "\"" . join("\n", @{$feat_ref->{$featid}->{'accessions'}->{'RAST'}}) . "\""
	    if ($feat_ref->{$featid}->{'accessions'}->{'RAST'});

	my @IMG_accs;
	push @IMG_accs, @{$feat_ref->{$featid}->{'accessions'}->{'IMG'}} if (defined $feat_ref->{$featid}->{'accessions'}->{'IMG'});
	push @IMG_accs, @{$feat_ref->{$featid}->{'accessions'}->{'JGI'}} if (defined $feat_ref->{$featid}->{'accessions'}->{'JGI'});
	my $IMG_id = "\"" . join("\n", @IMG_accs) . "\"";

	my $locus_tag = "\"" . join("\n", @{$feat_ref->{$featid}->{'accessions'}->{'locus_tag'}}) . "\""
	    if ($feat_ref->{$featid}->{'accessions'}->{'locus_tag'});



	# Annotation
	my %EC;
	my $RAST_ann = {};
	my $IMG_ann = {};
	my $cofactor = {};
	my $RAST_product;
	my $IMG_product;

	$RAST_ann = $feat_ref->{$featid}->{'annotation'}->{9}->{'RAST'} if (defined $feat_ref->{$featid}->{'annotation'}->{9}->{'RAST'});
	if (defined $RAST_ann->{'product'}) {
	    $RAST_product = join("\n", @{$RAST_ann->{'product'}}); }
	if (defined $RAST_ann->{'EC_number'}) { 
	    foreach my $ecr (@{$RAST_ann->{'EC_number'}}) {
		foreach my $ec (split(/\s+/, $ecr)) {
		    $EC{$ec} = 1; } } }
	$IMG_ann = $feat_ref->{$featid}->{'annotation'}->{9}->{'IMG'} if (defined $feat_ref->{$featid}->{'annotation'}->{9}->{'IMG'});
	if (defined $IMG_ann->{'product'}) {
	    $IMG_product = join("\n", @{$IMG_ann->{'product'}}); }
	if (defined $IMG_ann->{'EC_number'}) { 
	    foreach my $ecr (@{$IMG_ann->{'EC_number'}}) {
		foreach my $ec (split(/\s+/, $ecr)) {
		    $EC{$ec} = 1; } } }
	&get_cofactor_by_EC($dbh, \%EC, $cofactor);
	
	
	my $ev_ref = &get_evidence_for_feature($dbh, $featid);
	my %TIGR;
	my %TIGRdesc;
	my %PF;
	my %PFdesc;
	my $TMH;
	my $SP;
	my $OMP;
	my $mainrole;
	my $subrole;
	my %KO;
	my $KO_desc;
	my $CAZY = {};
	my $MEROPS = {};
	my $TC = {};
	my $COG;
	my $COG_desc;
	my $role = {};
	my $subsystem;
	my $COG_cat;
	my $SEED_subsys;

	foreach my $evrow(@$ev_ref) {
	    # 0 feat_id
	    # 1 feat_min
	    # 2 feat_max
	    # 3 program
	    # 4 ev_type
	    # 5 ev_accession
	    # 6 ev_min
	    # 7 ev_max
	    # 8 ev_length
	    # 9 score
	    if ($evrow->[4] eq "HMM") {
		my $hmm_ann = &get_hmm_annotation($dbh, $evrow->[5]);
		if ($hmm_ann->{$evrow->[5]}->{'ec_num'}) {
		    foreach my $ec (split(/\s+/,$hmm_ann->{$evrow->[5]}->{'ec_num'})) {
			$EC{$ec} = 1 } }
		if ($evrow->[5] =~ /^TIGR/) {
		    $TIGR{$evrow->[5]} = 1;
		    $TIGRdesc{$hmm_ann->{$evrow->[5]}->{'product'}} = 1;
		    &get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor, $role);
		} elsif ($evrow->[5] =~ /^PF/) {
		    $evrow->[5] =~ s/\.\d+//;
		    $PFdesc{$hmm_ann->{$evrow->[5]}->{'product'}} = 1;
		    &get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor, $role);
		    $PF{$evrow->[5]} = 1;
		} else {
		}
	    } elsif ($evrow->[4] eq "KO") {
		&get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor);
		$KO_desc = &get_KO_desc($dbh, $evrow->[5]);
		$KO{$evrow->[5]} = 1;
	    } elsif ($evrow->[4] eq "FIGFAM") {
		$SEED_subsys .= &get_SEED_subsystem($dbh, $evrow->[5]);
	    } elsif ($evrow->[4] eq "COG") {
		($COG, $COG_desc, $COG_cat) = &get_COG_and_cat($dbh, $evrow->[5]);
	    } elsif ($evrow->[4] eq "SP") {
		$SP .= "$evrow->[5]:";
	    } elsif ($evrow->[4] eq "TMH") {
		$TMH++;
	    } elsif ($evrow->[4] eq "CAZy") {
		&get_CAZY_from_CAZY($evrow->[5], $CAZY);
	    } elsif ($evrow->[4] eq "OMP") {
		$OMP .= "$evrow->[5]:";
	    }
	}

	# Now print it all out
	$output .= "$featid";
	$output .= "\t$RAST_id";
	$output .= "\t$PNNL_id";
	$output .= "\t$IMG_id";
	$output .= "\t$feat_type";
	$output .= "\t$locus_tag";
	$output .= "\t$contig";
	$output .= "\t$start";
	$output .= "\t$stop";
	$output .= "\t$strand";
	$output .= "\t$aa";
	$output .= "\t$subsystem";
	$output .= "\t\"" . join("\n", (keys %{$role->{'main'}})) . "\"";
	$output .= "\t\"" . join("\n", (keys %{$role->{'sub'}})) . "\"";
	$output .= "\t$COG_cat";
	$output .= "\t$SEED_subsys";
	$output .= "\t\"" . join(" ", (keys %$TC)) . "\"";
	$output .= "\t\"" . join(" ", (keys %KO)) . "\"";
	$output .= "\t$KO_desc";
	$output .= "\t\"$RAST_product\"";
	$output .= "\t\"$IMG_product\"";
	$output .= "\t\"" . join(" ", (keys %EC)) . "\"";
	$output .= "\t$COG";
	$output .= "\t$COG_desc";
	$output .= "\t\"" . join(" ", (keys %PF)) . "\"";
	$output .= "\t\"" . join("\n", (keys %PFdesc)) . "\"";
	$output .= "\t\"" . join(" ", (keys %TIGR)) . "\"";
	$output .= "\t\"" . join("\n", (keys %TIGRdesc)) . "\"";
	$output .= "\t$SP";
	$output .= "\t$TMH";
	$output .= "\t\"" . join("\n", (keys %{$cofactor->{'ec'}})) . "\"";
	$output .= "\t\"" . join("\n", (keys %{$cofactor->{'ev'}})) . "\"";
	$output .= "\t\"" . join("\n", (keys %$CAZY)) . "\"";
	$output .= "\t\"" . join("\n", (keys %$MEROPS)) . "\"";
	$output .= "\n";
    }
}

my $OUT;
if ($outfile){
    open $OUT, ">$outfile" || die "Can't open $outfile for write:$!";
} else {
    $OUT = \*STDOUT;
}

print $OUT &header(-type => 'text/plain',
		   -Content_disposition => "attachment; filename=$setname.tsv")
    . "Internal id\tRAST acc\tPNNL acc\tIMG acc\tFeat_type\tLocus_tag\tContig\tStart\tStop\tStrand\tAA\tSubsystem\tMainrole\tSubrole\tCOG Cat\tRAST subsystem\tTC\tKO\tKO description\tRAST product descriptor\tIMG product descriptor\tEC\tCOG\tCOG description\tPfam\tPfam description\tTIGRfam\tTIGRfam description\tSP\tTMH\tEC_Cofactor\tev_Cofactor\tCAZy\tMEROPS\n"
    . $output;

sub get_CAZY_from_CAZY {
    my $ev_acc = shift;
    my $CAZY = shift;
    
    $ev_acc =~ s/\.hmm$//;
    $CAZY->{$ev_acc} = 1;
}

sub get_MEROPS_cofactor_TC {
    my $dbh = shift;
    my $ev_acc = shift;
    my ($CAZY, $MEROPS, $TC, $cofactor, $role) = @_;

    my $MEROPSq = "select clan, family, subfamily from egad.ev_merops_link"
	. " WHERE ev_acc = \"$ev_acc\"";
    my $meropsr = $dbh->selectrow_arrayref($MEROPSq);
    if (defined $meropsr) {
	if ($meropsr->[2]) { 
	    $MEROPS->{$meropsr->[2]} = 1;
	} elsif ($meropsr->[1]) {
	    $MEROPS->{$meropsr->[1]} = 1;
	} elsif ($meropsr->[0]) {
	    $MEROPS->{$meropsr->[0]} = 1;
	}
#	$MEROPS->{join(':', @$meropsr)} = 1;
    }
    my $TCq = "select tc from egad.ev_tc_link"
	. " WHERE ev_acc = \"$ev_acc\"";
    my $tcr = $dbh->selectcol_arrayref($TCq);
    foreach my $r(@$tcr) { $TC->{$r} = 1 }

    my $cofq = "SELECT cofactor from egad.ev_cofactor_link where ev_acc=\"$ev_acc\"";
    my $cofr = $dbh->selectcol_arrayref($cofq);
    foreach my $r(@$cofr) { $cofactor->{'ev'}->{$r} = 1 }

    my $roleq = "SELECT m.name, s.name from egad.ev_role_link l, egad.subrole s, egad.mainrole m"
	. " WHERE l.ev_acc=\"$ev_acc\""
	. " AND s.id=l.subrole"
	. " AND m.id = left(l.subrole,2)";
    my $roler = $dbh->selectall_arrayref($roleq);
    foreach my $row(@$roler) {
	$role->{'main'}->{$row->[0]} = 1;
	$role->{'sub'}->{$row->[1]} = 1;
    }
}

sub get_hmm_annotation {
    my $dbh = shift;
    my $hmm_acc = shift;
    my $q = "SELECT hmm_acc, hmm_com_name AS product, ec_num, ko, localization, gene_sym, cazy, merops, tc_num"
	. " FROM egad.hmm2"
	. " WHERE hmm_acc = \"$hmm_acc\"";
    my $r = $dbh->selectall_hashref($q, 'hmm_acc');
    return $r
}

sub get_cofactor_by_EC {
    $dbh = shift;
    my $EC = shift;
    my $cofactor = shift;
    
    my $cofq = "SELECT cofactor from egad.ec_cofactor_link"
	. " WHERE ec=?";
    my $sth = $dbh->prepare($cofq);
    foreach my $ec (keys %$EC) {
	$sth->execute($ec);
	while (my $r = $sth->fetchrow_array) {
	    $cofactor->{'ec'}->{$r} = 1;
	}
    }
}
	    
sub get_SEED_subsystem {
    my ($dbh, $acc) = @_;
    my $q = "select category, subcategory, subsystem from egad.SEED s, egad.FIGFAM f where f.accession=\"$acc\" and s.function_id=f.description_id";
    my $r = $dbh->selectall_arrayref($q);
    my $string;
    foreach my $row(@$r) {
	$string .= join(":", @$row) . ";";
    }
    return $string;
}

sub get_KO_desc {
    my ($dbh, $ko) = @_;
    my $q = "select product, gene_sym, ec from egad.ko where ko=\"$ko\"";
    my $r = $dbh->selectall_arrayref($q);
    my $string;
    foreach my $row(@$r) {
	$string = $row->[0];
	$string .= " ($row->[1])" if ($row->[1]);
	$string .= " [$row->[2]]" if ($row->[2]);
    }
    return $string;
}

sub get_COG_and_cat {
    my ($dbh, $cog) = @_;
    my $q = "select description, category from egad.COG where accession=\"$cog\"";
    my $r = $dbh->selectall_arrayref($q);
    return $cog, $r->[0]->[0], $r->[0]->[1];
}
