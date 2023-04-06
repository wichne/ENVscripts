#!/usr/bin/perl

use strict;
use lib $ENV{ENVSCRIPTS};
use ENV;
use DBI;
use Getopt::Std;
use Bio::SeqIO;

my $opt = {};
&getopts('D:i:n:o:', $opt);
my $db = $opt->{D};
my $setid = $opt->{i};
my $setname = $opt->{n};
my $outfile = $opt->{o};

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
my $dbh = connect($opt);


my $OUT;
if ($outfile) {
    open $OUT, ">$outfile" || die "Can't open $outfile for write:$!";
} else {
    $OUT = \*STDOUT;
}

#my $setref = &get_seq_sets($dbh);
if ($setname && !$setid) { $setid = &set_name_to_id($dbh, $setname); }

my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS");

print $OUT join("\t", 'Set', 'ID', 'locus_tag', 'Type', 'SeqID', 'Min coord', 'Max coord', 'Strand', 'Phase', 'Product description', 'Gene', 'EC number', 'TCDB', 'CAZy', 'MEROPS', 'cofactor', 'KO', 'Mainrole', 'Subrole', 'Subsystem', 'COG', 'COG cat', 'SEED subsystem') . "\n";
foreach my $fid (sort {
    $protref->{$a}->{'location'}->{'seq_id'} <=> $protref->{$b}->{'location'}->{'seq_id'} ||
	$protref->{$a}->{'location'}->{'feat_min'} <=> $protref->{$b}->{'location'}->{'feat_min'}
		 }
		 keys %$protref) {

    # The sequence accession
    my $saccref = &get_accession_by_seq_id($dbh, $protref->{$fid}->{'location'}->{'seq_id'});
    my $sacc = $saccref->[0];
    my $locus_tag;
    # The location stuff
    my ($min, $max) = ($protref->{$fid}->{'location'}->{'feat_min'},
		       $protref->{$fid}->{'location'}->{'feat_max'});
    my $strand = $protref->{$fid}->{'location'}->{'strand'};
    if ($strand eq "+" || $strand == 1) { $strand = "+" }
    elsif ($strand eq "-" || $strand == -1) { $strand = "-" }
    elsif ($strand eq "0") { $strand = "." }
    else { $strand = "?" }

    my $phase = $protref->{$fid}->{'location'}->{'phase'};
    if (! defined $phase || $phase eq "") { warn "No phase reported for CDS $fid: defaulting to phase=0\n";
			    $phase = 0; }
    elsif ($phase !~ /^[012]$/) { warn "Invalid phase ('$phase') for CDS $fid: defaulting to phase=0\n";
				$phase = 0; }

    # The description stuff
    my $ID = "$db|$fid";
    my @alias;
    foreach my $src (keys %{$protref->{$fid}->{'accessions'}}) {
	my $acc = $protref->{$fid}->{'accessions'}->{$src};
	$locus_tag = $acc if ($src eq 'locus_tag');
    }
    my %ann;
    while (my ($qual,$vref) = each %{$protref->{$fid}->{'annotation'}}) {
	my ($k) = sort {$a<=>$b} keys %$vref;
	$ann{$qual} = $vref->{$k}->[0]->{value};
    }
    my $ev_ref = &get_evidence_for_feature($dbh, $fid);
    my $TMH;
    my $SP;
    my $OMP;
    my %KO;
    my $CAZY = {};
    my $MEROPS = {};
    my $TC = {};
    my $COG;
    my $COG_cat;
    my $COG_desc;
    my $role = {};
    my $subsystem;
    my $SEED_subsys;
    my $cofactor;
    foreach my $evrow(@$ev_ref) {
	if ($evrow->[4] eq "HMM") {
	    if ($evrow->[5] =~ /^TIGR/) {
		    &get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor, $role);
		} elsif ($evrow->[5] =~ /^PF/) {
		    &get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor, $role);
		} else {
		}
	} elsif ($evrow->[4] eq "KO") {
	    &get_MEROPS_cofactor_TC($dbh, $evrow->[5], $CAZY, $MEROPS, $TC, $cofactor);
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

    print $OUT join("\t", $setname, $fid, $locus_tag, $protref->{$fid}->{feat_type}, $sacc, $min, $max, $strand, $phase, $ann{'product'}, $ann{'gene'}, $ann{'ec'}, join(" ", keys %$TC), join(" ", keys %$CAZY), join(" ", keys %$MEROPS), join(" ", keys %{$cofactor->{'ev'}}), join(" ", keys %KO), join("::", keys %{$role->{'main'}}), join("::", keys %{$role->{'sub'}}), $subsystem, $COG, $COG_cat, $SEED_subsys );
    print $OUT "\n";

}


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

sub get_COG_and_cat {
    my ($dbh, $cog) = @_;
    my $q = "select description, category from egad.COG where accession=\"$cog\"";
    my $r = $dbh->selectall_arrayref($q);
    return $cog, $r->[0]->[0], $r->[0]->[1];
}
