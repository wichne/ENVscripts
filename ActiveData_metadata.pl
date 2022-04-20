#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

&getopts('f:');
our $dbh = connect({'u' => 'access',
		   'p' => 'access',
		   'D' => 'hotlake_ucc',});

our $MAINROLE = &get_mainroles($dbh);
our $SUBROLE = &get_subroles($dbh);

if ($opt_f) {
    open (IN, $opt_f) || die "Can't open $opt_f for reading: $!\n";
    my @accs;
    while (my $l = <IN>) {
	chomp $l;
	my @f = split/[\s\,]+/, $l;
	push @accs, $f[0];
    }
    my $featref = get_feature_id_by_accession($dbh, @accs);
    my @fids;
    foreach my $ref(values %$featref) { push @fids, $ref->{'feature_id'} }
    my $protref = get_features_by_feature_id($dbh, @fids);
    &process_protref($protref);
} else {
    $set_ref = get_set_names($dbh);
    while (my($id, $ref) = each %$set_ref) {
	push @set_ids, $id;
	push @set_names, $ref->{'name'};
    }
    
    for (my $i=0; $i<@set_ids; $i++) {
	my $set_name=$set_names[$i];
	my $setid = $set_ids[$i];
	my $protref = &get_seq_features_by_set_id($dbh, $setid, "CDS", "tRNA", "rRNA", "ncRNA");
	
	&process_protref($protref);
    }
}

sub process_protref {
    my ($protref, $fid) = @_;

    foreach my $fid (keys %$protref) {
	my $feat_type = $protref->{$fid}->{'feat_type'};
	my $locus_tag = $protref->{$fid}->{'accessions'}->{'locus_tag'};
	# get KO
	my $KO_ref = get_evidence_for_feature($dbh, $fid, "KO");
	my $KO = $KO_ref->[0]->[5];

	# get role from best evidence
	my $q = "SELECT distinct subrole, LEFT(subrole,2) AS mainrole, fe.ev_accession, fe.score"
	    . " FROM egad.ev_role_link erl, feature_evidence fe"
	    . " WHERE fe.feature_id=$fid"
	    . " AND erl.ev_acc=fe.ev_accession"
	    . " ORDER BY fe.score DESC";
	my $role_ref = $dbh->selectall_arrayref($q);
	my $mainrole = $MAINROLE->{$role_ref->[0]->[1]}->{'name'} ? $MAINROLE->{$role_ref->[0]->[1]}->{name} : "";
	my $subrole = $SUBROLE->{$role_ref->[0]->[0]}->{'name'} ?  $mainrole . ":" . $SUBROLE->{$role_ref->[0]->[0]}->{name} : "";

	print "$fid|$locus_tag\t$feat_type\t$KO\t$mainrole\t$subrole\n";
    }
}    
