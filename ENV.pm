#!/usr/bin/perl
package ENV;
require Exporter;
use Carp;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use strict;
use vars qw (@ISA @EXPORT);
use DBI;
@ISA = qw(Exporter);
@EXPORT = qw(connect get_sequences load_SeqFeature load_Sequence load_feature_annotations delete_feature_annotations get_sequence_by_feature_id get_accession_by_seq_id get_feature_ids_by_seq_id get_features_by_feature_id get_sequence_by_seq_id get_seq_id_by_seq_accession link_seq_to_set get_seq_sets load_sequence_sets load_PrimarySeq load_sequence_SeqFeature load_sequences load_sequence_accessions load_sequence_annotations get_sequence_by_accession set_name_to_id get_seq_features_by_set_id get_sequences_by_set_id get_feature_id_by_accession get_mainroles get_subroles set_id_to_seq_ids seq_id_to_SeqObj add_features_to_SeqObj get_evidence_for_feature get_ev_role_by_feature_id get_accession_by_feature_id get_features_by_seq_id load_seq_feat_mappings set_id_to_set_name get_locus_tag_by_feature_id get_set_names seq_id_to_set_name is_current update_product insert_feature_accessions update_feature_mapping set_name_to_set_id);

sub connect {
    my $argref = shift;
    my $dbh;
    my ($user, $db, $password);
    if (defined $argref->{'u'} && $argref->{'u'}) { $user = $argref->{'u'} }
    elsif (defined $argref->{'U'} && $argref->{'U'}) { $user = $argref->{'U'} }
    elsif (defined $argref->{'user'} && $argref->{'user'}) { $user = $argref->{'user'} }
    else { $user = $ENV{USER}; }
    
    if (defined $argref->{'p'}) { $password = $argref->{'p'} }
    elsif (defined $argref->{'P'}) { $password = $argref->{'P'} }
    elsif (defined $argref->{'pswd'}) { $password = $argref->{'pswd'} }
    
    if (defined $argref->{'D'}) { $db = $argref->{'D'} }

    my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';
    $dbh = DBI->connect("dbi:mysql:host=$host;database=$db", $user, $password);
    if (! $dbh) { die "Couldn't connect with $host / $db / $user / $password\n";}
    print STDERR "Connected to $host:$db as $user\n";
    return $dbh;

}

sub get_set_names {
    my $dbh = shift;
    my $set_q = "select set_id, name from sequence_sets where is_current=1";
    return $dbh->selectall_hashref($set_q, 'set_id');
}

sub get_all_sequences {
    my $dbh = shift;
    my $limit_ref = shift;

    my $sequence_q = "SELECT s.seq_id, sequence, seq_accession"
	. " FROM sequences s, sequence_accessions a"
	. " WHERE s.seq_id = a.seq_id"
	. " AND a.seq_acc_description = 'accession_number'";
#    my $sth = $dbh->prepare($sequence_q);
#    $sth->exectute;
#    while ($row = $sth->
    my $r = $dbh->selectall_hashref($sequence_q, 'seq_id');
    return $r;
}

sub get_sequence_by_accession {
    my $dbh = shift;
    my $accession = shift;

    my $sequence_q = "SELECT s.seq_id, sequence, seq_accession"
	. " FROM sequences s, sequence_accessions a"
	. " WHERE s.seq_id = a.seq_id"
	. " AND a.seq_accession = \"$accession\"";
    my $r = $dbh->selectall_hashref($sequence_q, 'seq_id');
    if (scalar(keys %$r) > 1) { die "ERROR: Sequence accession '$accession' has multiple seq_ids in database. Dying...\n";}
    foreach my $sid (keys %$r) {
	return Bio::Seq->new(-seq => $r->{$sid}->{'sequence'},
			     -id => $sid,
			     -accession => $r->{$sid}->{'seq_accession'});
    }
}

sub get_sequence_by_seq_id {
    my $dbh = shift;
    my @seq_ids = @_;
    
    my $data = {};
    my $sequence_q = "SELECT s.seq_id, sequence"
	. " FROM sequences s"
	. " WHERE s.seq_id in (". join (",",@seq_ids) . ")";
    my $r = $dbh->selectall_hashref($sequence_q, 'seq_id') or warn "sequence_q $DBI::errstr";
    foreach my $sid (keys %$r) {
	$data->{$sid} = Bio::Seq->new(-seq => $r->{$sid}->{'sequence'},
				      -id => $sid,
				      -accession => $r->{$sid}->{'seq_accession'});
	# grab annotation for header
	my $annotation_q = "SELECT sa.data_type, sa.value"
	    . " FROM sequence_annotations sa" #, env.data_types d"
	    . " WHERE sa.seq_id=$sid";
#	. " AND sa.data_type_id=d.id";
	my $annotations = $dbh->selectall_arrayref($annotation_q) or warn "$annotation_q $DBI::errstr";
	my $desc;
	foreach (@$annotations) {
	    $desc .= "$_->[0]=$_->[1];";
	}
	$data->{$sid}->desc($desc);
    }
    return $data;
}

sub get_seq_id_by_seq_accession {
    my $dbh = shift;
    my $accession = shift;

    my $seq_id_q = "SELECT distinct x.seq_id FROM sequence_accessions x, sequences s"
	. " WHERE seq_accession = \"$accession\""
	. " AND s.seq_id=x.seq_id"
	. " AND s.iscurrent=1";

    my $r = $dbh->selectall_arrayref($seq_id_q);
    return $r->[0]->[0];
}

sub get_accession_by_seq_id {
    my $dbh = shift;
    my $seq_id = shift;

    my $seq_id_q = "SELECT seq_accession FROM sequence_accessions"
	. " WHERE seq_id = \"$seq_id\"";

    my $r = $dbh->selectcol_arrayref($seq_id_q);
    return $r;
}

sub get_accession_by_feature_id {
    my $dbh = shift;
    my $f_id = shift;
    my $src = shift;

    my $f_id_q = "SELECT accession FROM feature_accessions"
	. " WHERE feature_id = \"$f_id\"";
    if ($src) {
	$f_id_q .= " AND source=\"$src\"";

	my $r = $dbh->selectcol_arrayref($f_id_q);
	return $r;
    }
}

sub get_locus_tag_by_feature_id {
    my $dbh = shift;
    my $f_id = shift;

    my $f_id_q = "SELECT accession FROM feature_accessions"
	. " WHERE feature_id = \"$f_id\""
	. " AND source = \"locus_tag\"";

    my $r = $dbh->selectcol_arrayref($f_id_q);
    return $r;
}

sub get_feature_id_by_accession {
    my $dbh = shift;
    my $accession = shift;
    if ($accession =~ /^\s*$/) { warn "Empty accession requested ($accession)\n" }
    my @a = split(/\|/, $accession);
    my $f_id_q = "SELECT distinct feature_id FROM feature_accessions"
	. " WHERE accession in (\"" . join("\",\"", @a) . "\")";
    my $r = $dbh->selectcol_arrayref($f_id_q);
    if (defined $r) {
	if (@$r > 1) {
	    warn "$accession yielded more than one feature_id: " . join(" ", @$r) 
		. ". Returning the first.\n"; }
	return $r->[0];
    } else { return }
}

sub load_SeqFeature {
    my $dbh = shift;
    my $seq_id = shift;
    my $feato = shift; # Bio::SeqFeature::Generic object
    my $seqobj = shift; # as above, except for the parent DNA sequence
    my $SO_term = shift;
    my $source = shift;
    my $prefix = shift;

    if (!defined $seqobj) {
	my $seqstruct = get_sequence_by_seq_id($dbh, $seq_id);
	$seqobj = Bio::Seq->new(-seq => $seqstruct->{$seq_id}->{sequence},
				-alphabet => 'dna');
    }

    # could check for split location with
    # if ( $feato->location->isa('Bio::Location::SplitLocationI')
    # see BioPerl or CPAN for more info
    # start just gives the start. If we want to see the extendable, use to_FTstring or start_pos_type
    my $start = $feato->start;
    my $min_partial = $feato->location->start_pos_type eq "BEFORE" ||
	$feato->location->min_start != $feato->location->max_start ? 1 : 0; 
    
    my $end = $feato->end;
    my $max_partial = $feato->location->end_pos_type eq "AFTER" ||
	$feato->location->min_end != $feato->location->max_end ? 1 : 0; 
    my $pseudo = $feato->has_tag('pseudo') ? 5 : 0;

    # since peseudogenes don't have a CDS entry...
    my $feat_type = $feato->primary_tag eq "gene" && $pseudo > 0 ? "CDS" : $feato->primary_tag;

    my $prot = "";
    my $phase = 0;
    my @all_tags = $feato->get_all_tags;
    if ($feato->primary_tag eq "CDS") {
#    if ($feat_type_id == 1) { # if this is a translated gene...
	if (grep /\btranslation\b/, @all_tags) {
	    ($prot) = $feato->get_tag_values("translation");
	} else {
	    my $protobj = protein_from_coords($seqobj, $feato);
	    $prot = $protobj->seq;
	}
	if (defined $feato->frame) {
	    $phase = $feato->frame;
	} else {
	    my @codon_start = $feato->annotation->get_Annotations('codon_start');
	    if (defined $codon_start[0]) {
		$phase = $codon_start[0]->value - 1;
	    }
	}
    }
    
    $prot =~ s/\*$//;
    my $sequence_features_i = "INSERT sequence_features"
	. " (feat_type_id, feat_type, product, inserted_by, date_inserted)"
	. " SELECT id, \"$feat_type\", \"$prot\", USER(), NOW()"
	. " FROM INSDC.feature_key"
	. " WHERE feature_key = \"$feat_type\"";
#    print $sequence_features_i, "\n";
    my $row1 = $dbh->do($sequence_features_i);
    if (! defined $row1) {
	warn "\n!!!Couldn't execute '$sequence_features_i': $DBI::errstr\n";
	return;
    }
    
    my $strand = $feato->strand;
    my $translation_coords = $feato->location->to_FTstring;
    my $feat_id = $dbh->last_insert_id("%", "%", "", "");
    my $seq_feat_mappings_i = sprintf "INSERT seq_feat_mappings"
	. " (seq_id, feature_id, feat_min, feat_max, strand,"
	. " phase, min_partial, max_partial, pseudo, translation_coords)"
	. " VALUES ($seq_id, $feat_id, $start, $end, '$strand',"
	. " '$phase', $min_partial, $max_partial, $pseudo, '$translation_coords')";
#    print $seq_feat_mappings_i, "\n";
    my $row2 = $dbh->do($seq_feat_mappings_i);
    if (! defined $row2) {
	warn "\n!!!Couldn't execute '$seq_feat_mappings_i'\n";
	return;
    }
    
    # insert an accession. Use the locus_tag, if present
    my $accession = $feato->primary_id;
    if ($feato->has_tag('locus_tag')) {
	$accession = [$feato->get_tag_values('locus_tag')]->[0];
    }
    elsif ($feato->display_name) {
	$accession = $feato->display_name;
    }

    my $feat_acc_i = "INSERT feature_accessions"
	. " (feature_id, source, prefix, accession)"
	. " VALUES ($feat_id, \"$source\", \"$prefix\", \"$accession\")" ;
    my $row3 = $dbh->do($feat_acc_i);
    
    # insert annotation
    if (grep /^product$/, @all_tags) {
	my $feat_ann_i = "INSERT into feature_annotations"
	    . " (feature_id, data_type_id, value, source, rank)"
	    . " VALUES ($feat_id, 66, \"" . [$feato->get_tag_values('product')]->[0] . "\", \"$source\", 10)";
	my $row4 = $dbh->do($feat_ann_i);
    }
    
    return $feat_id;
}

sub protein_from_coords {
    my ($seqobj, $featobj) = @_;
    if (ref $seqobj eq "SCALAR") {
	# we've been passed a seq_id, so get the seqobj
#	$seqobj = get_sequence_by_seq_id($dbh, $seqobj);
    }

    # Cut gene sequence from contig
    my $featseqobj = $seqobj->trunc($featobj->location);
    
    # Translate sequences to aa if applicable
    my $complete = $featobj->location->start_pos_type ne "EXACT" ||
	$featobj->location->end_pos_type ne "EXACT" ? 0 : 1; 
    my ($table) = $featobj->get_tag_values("transl_table");
    $table = 11 if (! $table);
    my $protobj = $featseqobj->translate(-complete => $complete,
					 #-frame    => $phase,
					 -codontable_id => $table);
    return ($protobj);
}

sub delete_feature_annotations {
    my $dbh = shift;
    my $feat_id = shift;
    my $cond_ref = shift;
    my $feature_annotations_d = "DELETE FROM feature_annotations"
	. " WHERE feature_id=$feat_id";
    if (defined %$cond_ref) {
	foreach my $cond(keys %$cond_ref) {
	    $feature_annotations_d .= " AND $cond = \'$cond_ref->{$cond}\'";
	}
    }
    my $d = $dbh->do($feature_annotations_d);
    if (! defined $d) {
	warn "Couldn't execute '$feature_annotations_d'";
	return;
    }
}
sub load_feature_annotations {
    my $dbh = shift;
    my $feat_id = shift;
    my $annObj = shift;
    my $source = shift;
    my $rank = shift;

    my $feature_annotation_i = "INSERT feature_annotations"
	. " (feature_id, data_type_id, value, edit_by, date_edit, rank, source)"
	. " SELECT $feat_id, d.id, ?, USER(), NOW(), ?, ?"
	. " FROM INSDC.qualifier d"
	. " WHERE d.qualifier=?";

    foreach my $data_type ($annObj->get_all_annotation_keys()) {
	my @values = $annObj->get_Annotations($data_type);
	for (my $i = 0; $i < @values; $i++) {
	    if ($values[$i]->value eq "") { next }
	    if ($data_type eq "db_xref") {
		my ($type, $value) = split/\:/, $values[$i]->value;
		$dbh->do($feature_annotation_i, {}, ($value, 0, $type, $source)) or warn $DBI::errstr;
	    } else {
		if (! $rank) { $rank = 9 }
		$dbh->do($feature_annotation_i, {}, ($values[$i]->value, $rank, $source, $data_type)) or warn $DBI::errstr;
	    }
	}
    }	    
}

sub load_PrimarySeq {
    my $dbh = shift;
    my $seqo = shift;

    my $fail = 0;
    # error checking
    if (! defined $dbh ) { warn "No dbh provided"; return }
    if (! defined $seqo) { warn "No seqo provided"; return }
    if (ref $seqo !~ /HASH/) { warn "seqo is not a hash reference! ($seqo)"; return }

    my $seq = $seqo->seq;
    if (! $seq) { warn "No sequence in object"; return }
    my $seq_id = &load_sequences($dbh, $seq);

    my $acc = $seqo->primary_id ? $seqo->primary_id : 
	$seqo->display_name ? $seqo->display_name : undef;
    if (! defined $acc) { warn "No accession found in object"; return }
    print STDERR "$acc loaded as seq_id $seq_id\n";
    my ($prefix, $source);
    my $acc_ref = {};
    my @acc = split/\|/, $acc;
    foreach my $a(@acc) {
	$source = "undetermined";
	if ($a eq "quiver") { next }
	if ($a =~ /^\d{13}$/) { $source = "uid" }
	elsif ($a =~ /^[A-Z]{4}\w\d{2}T[FR]/) { $source = "seq_name" }
	$acc_ref->{$source} = $acc;
    }
    my $row2 = &load_sequence_accessions($dbh, $seq_id, {$source => $acc});
    if ($row2 > 0) {
	$fail = 1;
    }

    my $desc = $seqo->desc;
    my %DESC; # holds key/value pairs inside brackets
    # NCBI: identify organism and remove from description
    if ($desc =~ /(\{\s*(.+);.*\}\s*)/) {
	$DESC{'organism'} = $2;
	$desc =~ s/\Q$1\E//;
    }
    # identify other key/value pairs and remove from description
    # Look for Genbank format:
    if ($desc =~ /\[\w+=[^\]]+\]/) {
	while ($desc =~ /(\[(\w+)=([^\]]+)\])/g) {
	    $DESC{$2}=$3;
	    $desc =~ s/\Q$1\E\s*//;
	}
    }
    # look for JTC format
    elsif ($desc =~ /\/\w+=\S+/) {
	while ($desc =~ /(\/(\w+)=(\S+))/) {
	    $DESC{$2} = $3;
	    $desc =~ s/\Q$1\E\s*//;
	}
    }
    # what's left is description
    $DESC{'description'} = $desc if ($desc =~ /\S+/);
    my $row4 = &load_sequence_annotations($dbh, $seq_id, \%DESC);
    if ($row4 > 0) {
	$fail = 1;
    }

    if ($fail) {
	warn "Rolling back insertion due to failure.";
	$dbh->rollback;
	return;
    } else {
#	$dbh->commit;
	return ($seq_id);
    }
}

sub load_sequence_SeqFeature {
    my $dbh = shift;
    my $seqo = shift; # Bio::SeqFeatureI object

    my $fail = 0;
    # error checking
    if (! defined $dbh ) { warn "No dbh provided"; return }
    if (! defined $seqo) { warn "No seqo provided"; return }
    if (ref $seqo !~ /HASH/) { warn "seqo is not a hash reference! ($seqo)"; return }

    my $seq = $seqo->seq->seq;
    if (! $seq) { warn "No sequence in object"; return }
    my $seq_id = &load_sequences($dbh, $seq);

    my $acc = $seqo->seq_id ? $seqo->seq_id :
	$seqo->primary_id ? $seqo->primary_id : 
	$seqo->display_name? $seqo->display_name : undef;
    if (! defined $acc) { warn "No accession found in object"; return }
    # we should always have a primary_id
    my $source = "undetermined";
    if ($acc =~ /^\d{13}$/) { 
	$source = "uid";
    }
    elsif ($acc =~ /^[A-Z]{4}\w\d{2}T[FR]/) { $source = "seq_name" }
    elsif (defined $seqo->source_tag) { $source = $seqo->source_tag }
    
    my $row2 = &load_sequence_accessions($dbh, $seq_id, {$source => $acc});
    if ($row2 > 0) {
	$fail = 1;
    }
    
    my $desc = $seqo->seq->description;
    my %DESC; # holds key/value pairs inside brackets
    # NCBI: identify organism and remove from description
    if ($desc =~ /(\{\s*(.+);.*\}\s*)/) {
	$DESC{'organism'} = $2;
	$desc =~ s/\Q$1\E//;
    }
    # identify other key/value pairs and remove from description
    # Look for Genbank format:
    if ($desc =~ /\[\w+=[^\]]+\]/) {
	while ($desc =~ /(\[(\w+)=([^\]]+)\])/g) {
	    $DESC{$2}=$3;
	    $desc =~ s/\Q$1\E\s*//;
	}
    }
    # look for JTC format
    elsif ($desc =~ /\/\w+=\S+/) {
	while ($desc =~ /(\/(\w+)=(\S+))/) {
	    $DESC{$2} = $3;
	    $desc =~ s/\Q$1\E\s*//;
	}
    }
    # what's left is description
    $DESC{'description'} = $desc if ($desc =~ /\S+/);

    my $row4 = &load_sequence_annotations($dbh, $seq_id, \%DESC);
    if ($row4 > 0) {
	$fail = 1;
    }

    if ($fail) {
	warn "Rolling back insertion due to failure.";
	$dbh->rollback;
	return;
    } else {
#	$dbh->commit;
	return ($seq_id);
    }
}

sub load_sequences {
    my $dbh = shift;
    my $seq = shift;
    my $len = length($seq);

    my $sequences_i = "INSERT sequences"
	. " (sequence, date_inserted, inserted_by, seq_length, iscurrent)"
	. " VALUES (\"$seq\", NOW(), USER(), $len, 1)";
    my $row1 = $dbh->do($sequences_i);
    if (! defined $row1) {
	warn "Couldn't execute '$sequences_i'";
	return;
    }
    my $seq_id = $dbh->{'mysql_insertid'}
    or die "no insert id?";
    return $seq_id;
}

sub get_sequence_by_feature_id {
    use Bio::Perl;
    my $dbh = shift;
    my @feature_ids = @_;

    # make a temporary table of all the seq_feat mappings
    my $temp_table_q = "CREATE TEMPORARY TABLE tmp_feat_maps"
	. " LIKE seq_feat_mappings";
    $dbh->do($temp_table_q) || warn "Couldn't create temporary table :: $DBI::errstr";

    # populate it with all the relevant data
    my $temp_table_i = "INSERT INTO tmp_feat_maps"
	. " SELECT * from seq_feat_mappings WHERE feature_id IN (" . join (",", @feature_ids) . ")";
    $dbh->do($temp_table_i) ||  warn "Population of temporary table failed : $DBI::errstr";

    # now grab all the sequences
    my $subseq_q = "SELECT t.*, substring(s.sequence, t.feat_min, (t.feat_max - t.feat_min + 1)) AS gene_seq"
	. " FROM tmp_feat_maps t, sequences s"
	. " WHERE s.seq_id = t.seq_id";
    my $sequences = $dbh->selectall_hashref($subseq_q, "feature_id");

    # if the strand is -1, we need to reverse complement
    while (my ($id, $ref) = each %$sequences) {
	if ($ref->{'strand'} == -1) {
	    $ref->{'gene_seq'} = revcom(uc($ref->{'gene_seq'}))->seq;
	}
    }

    # grab annotation for header
    my $annotation_q = "SELECT fa.feature_id, q.qualifier, fa.value"
	. " FROM feature_annotations fa, tmp_feat_maps t, INSDC.qualifier q"
	. " WHERE fa.feature_id=t.feature_id"
	. " AND q.id=fa.data_type_id";
    my $annotations = $dbh->selectall_arrayref($annotation_q);

    foreach (@$annotations) {
	$sequences->{$_->[0]}->{'annotation'}->{$_->[1]} = $_->[2];
    }

    return $sequences;
}

sub get_features_by_seq_id {
    my $dbh = shift;
    my $seq_id = shift;
    my $feat_id_ref = &get_feature_ids_by_seq_id($dbh, $seq_id);
    return (&get_features_by_feature_id($dbh, @$feat_id_ref));
}

sub get_features_by_feature_id {
    my $dbh = shift;
    my @feature_ids = @_;
    if (!@feature_ids) { warn "Hey! No feature_ids provided to ENV::get_features_by_feature_id!!\n"; return;}

    # grab all the sequences
    my $prot_q = "SELECT feature_id, feat_type, product"
	. " FROM sequence_features"
	. " WHERE feature_id in (" . join(",",@feature_ids) . ")"
	;

    my $sequences = $dbh->selectall_hashref($prot_q, "feature_id");
    my $mq = "select m.feature_id, l.set_id, m.seq_id, feat_min, feat_max, strand,"
	. " phase, min_partial, max_partial, pseudo, translation_coords"
	. " FROM seq_feat_mappings m, sequences s, seq_set_link l"
	. " WHERE m.feature_id in (" . join(",", @feature_ids) . ")"
	. " AND m.seq_id=l.seq_id"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1";
    my $msth = $dbh->prepare($mq);
    my $r1 = $msth->execute;
    if (! defined $r1) { die "Died on \n$mq\n" }
    while (my $rv = $msth->fetchrow_hashref()) {
	# Allow lookup of seq_id by feature_id and set_id:
	$sequences->{$rv->{feature_id}}->{'seq_id'}->{$rv->{set_id}} = $rv->{seq_id};
	$sequences->{$rv->{feature_id}}->{'location'}->{$rv->{seq_id}} = $rv;	
    }

    # grab annotation for header
    my $annotations = get_feature_annotations_by_feature_id($dbh, @feature_ids);
#! Right now this annotation structure is different than that returned by
#! get_seq_features_by_set_id.
    foreach my $a(@$annotations) {
	push @{$sequences->{$a->[0]}->{'annotation'}->{$a->[3]}->{$a->[4]}->{$a->[1]}}, $a->[2];
    }

    # grab the accessions
    my $accessions_q = "SELECT feature_id, source, prefix, accession"
	. " FROM feature_accessions"
	. " WHERE feature_id in (" . join (",",@feature_ids) . ")";
    my $accessions = $dbh->selectall_arrayref($accessions_q);

    foreach my $a(@$accessions) {
	my $acc = $a->[2] ? "$a->[2]|$a->[3]" : $a->[3];
	push @{$sequences->{$a->[0]}->{'accessions'}->{$a->[1]}}, $acc; 
    }

    # grab the set 
    my $set_q = "SELECT m.feature_id, s.set_id, s.name"
	. " FROM sequence_sets s, seq_set_link l, seq_feat_mappings m"
	. " WHERE m.feature_id in (" . join(",",@feature_ids) . ")"
	. " AND l.seq_id=m.seq_id"
	. " AND l.set_id=s.set_id"
	. " AND s.is_current=1";
    my $sets = $dbh->selectall_arrayref($set_q);
    foreach my $s(@$sets) {
	$sequences->{$s->[0]}->{'set'}->{'id'} = $s->[1];
	$sequences->{$s->[0]}->{'set'}->{'name'} = $s->[2];
    }

    return $sequences;

}

sub get_feature_annotations_by_feature_id {
    my $dbh = shift;
    my @feature_ids = @_;
    if (!@feature_ids) { warn "Hey! No feature_ids provided to ENV::get_features_by_feature_id!!\n"; return;}
    # grab annotation for header
    my $annotation_q = "SELECT fa.feature_id, d.qualifier, fa.value, fa.rank, fa.source"
	. " FROM feature_annotations fa, INSDC.qualifier d"
	. " WHERE feature_id in (" . join(",",@feature_ids) . ")"
	. " AND d.id=fa.data_type_id"
	. " ORDER BY rank";
    my $annotations = $dbh->selectall_arrayref($annotation_q);
    return $annotations;
}

sub get_ev_role_by_feature_id {
    my $dbh = shift;
    my $fid = shift;
    my $q = "SELECT distinct subrole, LEFT(subrole,2) AS mainrole FROM egad.ev_role_link erl, feature_evidence fe"
	. " WHERE fe.feature_id=$fid"
	. " AND erl.ev_acc=fe.ev_accession";
    my $r = $dbh->selectall_arrayref($q);
    my %R = ('mainrole' => [],
	     'subrole' => []);
    foreach my $row (@$r) {
	push @{$R{'subrole'}}, $row->[0];
	push @{$R{'mainrole'}}, $row->[1];
    }
    return \%R;    
}

sub get_feature_ids_by_seq_id {
    my $dbh = shift;
    my $seq_id = shift;
    
    my $feat_id_q = "SELECT m.feature_id"
	. " FROM seq_feat_mappings m, sequence_features f"
	. " WHERE seq_id = $seq_id"
	. " AND f.feature_id=m.feature_id AND is_current=1";
    return $dbh->selectcol_arrayref($feat_id_q);
}

# Populates seq_set_link to link seq_ids to a set_id
sub link_seq_to_set {
    my $dbh = shift;
    my $set_id = shift;
    my @seq_id = @_;
    foreach my $seq_id(@seq_id) {
	my $seq_set_link_i = "INSERT seq_set_link (seq_id, set_id)"
	    . " VALUES ($seq_id, $set_id)";
	my $row = $dbh->do($seq_set_link_i);
	if (! defined $row) {
	    warn "Couldn't execute '$seq_set_link_i': $dbh->errstr";
	} else {
	    print STDERR "seq_id $seq_id was linked to set_id $set_id\n";
	}
    }
}

# Loads the sequence_sets table
sub load_sequence_sets {
    my $dbh = shift;
    my $name = shift;
    my $desc = shift;
    my $ncbi_proj_id =shift;

    my $sequence_sets_i = "INSERT sequence_sets (name, description, ncbi_proj_id)"
	. " VALUES (\"$name\", \"$desc\", \"$ncbi_proj_id\")";
    my $row = $dbh->do($sequence_sets_i);
    if (! defined $row) {
	warn "Couldn't execute '$sequence_sets_i'";
    }
    my ($set_id) = $dbh->selectrow_array("select last_insert_id()") or die "no insert id?";
    return $set_id;
    
}

sub get_seq_sets {
    my $dbh = shift;
    my $set_q = "SELECT set_id, name, description"
	. " FROM sequence_sets WHERE is_current=1";
    my $r = $dbh->selectall_hashref($set_q, 'set_id');
    return $r;
}

sub load_sequence_accessions {
    my $dbh = shift;
    my $seq_id = shift;
    my $acc_ref = shift;
    my $fails = 0;
    while (my ($source, $acc) = each %$acc_ref) {
	my $sequence_accessions_i = "INSERT sequence_accessions"
	    . " (seq_id, seq_accession, seq_acc_source, seq_acc_description)"
	    . " VALUES ($seq_id, \"$acc\", \"$source\", 'accession')";
	my $r = $dbh->do($sequence_accessions_i);
	if (! defined $r) {
	    $fails += 1;
	    warn "Couldn't execute $sequence_accessions_i: $DBI::errstr";
	}
    }
    return $fails;
}

sub load_sequence_annotations {
    my $dbh = shift;
    my $seq_id = shift;
    my $annotations_ref = shift;
    my $fails = 0;
    while (my ($key, $value) = each %$annotations_ref) {
	if ($key !~ /.+/ || $value !~ /.+/) { next }
	my $sequence_annotations_i = "INSERT sequence_annotations (seq_id, data_type, value)"
	    . " VALUES ($seq_id, \"$key\", \"$value\")";
#	    . " SELECT $seq_id, id, \"$value\""
#	    . " FROM env.data_types"
#	    . " WHERE name='$key'";
	my $r = $dbh->do($sequence_annotations_i);
	if (! defined $r) { $fails++; warn "Couldn't execute $sequence_annotations_i: $DBI::errstr";}
    }
    return $fails;
}

sub get_seq_features_by_set_id {
    my $dbh = shift;
    my $set_id = shift;
    my @feat_type = @_;

    my $sq = "select f.feature_id, product, feat_type"
	. " FROM sequence_features f, seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND f.feature_id=m.feature_id AND f.is_current=1"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1"
	;
    if (@feat_type) {
	$sq .= " AND feat_type in ('" . join("', '", @feat_type) . "')";
    }
    my $r = $dbh->selectall_hashref($sq, 'feature_id');
    if (! defined $r) { return }
    
    my $mq = "select m.feature_id, m.seq_id, feat_min, feat_max, strand,"
	. " phase, min_partial, max_partial, pseudo, translation_coords"
	. " FROM seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND s.seq_id = m.seq_id"
	. " AND s.iscurrent=1"
	;
    my $rv = $dbh->selectall_hashref($mq, 'feature_id');
    foreach my $fid (keys %$rv) {
	if (defined $r->{$fid}) {
	    $r->{$fid}->{'location'} = $rv->{$fid};
	}
    }

    my $nq = "select a.feature_id, qualifier, value, rank, source"
	. " FROM feature_annotations a, seq_set_link l, seq_feat_mappings m,"
	. " INSDC.qualifier q, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND a.feature_id=m.feature_id"
	. " AND q.id=a.data_type_id"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1"
	;
    my $nsth = $dbh->prepare($nq);
    $nsth->execute();
    while (my @rv = $nsth->fetchrow_array) {
	if (!$rv[3]) { $rv[3] = 9;}
	if (defined $r->{$rv[0]}) {
	    push @{$r->{$rv[0]}->{'annotation'}->{$rv[1]}->{$rv[3]}}, {'value' => $rv[2],
								       'source' => $rv[4]};
	}
    }

    my $rq = "SELECT fe.feature_id, subrole, LEFT(subrole,2) AS mainrole"
	. " FROM egad.ev_role_link erl, feature_evidence fe, seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND fe.feature_id=m.feature_id"
	. " AND erl.ev_acc=fe.ev_accession"
	. " AND s.seq_id=m.seq_id and s.iscurrent=1"
	. " GROUP BY fe.feature_id, subrole, mainrole";
    my $rsth = $dbh->prepare($rq);
    $rsth->execute();
    while (my @rv = $rsth->fetchrow_array) {
	if (defined $r->{$rv[0]}) {
	    push @{$r->{$rv[0]}->{'role'}->{'mainrole'}}, $rv[2];
	    push @{$r->{$rv[0]}->{'role'}->{'subrole'}}, $rv[1];
	}
    }
    
    my $hcq = "SELECT e.feature_id, cofactor"
	. " FROM feature_evidence e, egad.ev_cofactor_link c, seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND e.feature_id=m.feature_id"
	. " AND (e.ev_accession=c.ev_acc OR LEFT(e.ev_accession,7)=c.ev_acc)"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1"
	. " GROUP by e.feature_id, cofactor";
    my $hcsth = $dbh->prepare($hcq);
    $hcsth->execute();
    while (my @rv = $hcsth->fetchrow_array) {
	if (defined $r->{$rv[0]}) {
	    push @{$r->{$rv[0]}->{'cofactor'}->{'hmm'}}, $rv[1];
	}
    }
    
    my $ecq = "SELECT a.feature_id, cofactor"
	. " FROM feature_annotations a, egad.ec_cofactor_link c, seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND a.feature_id=m.feature_id"
	. " AND a.data_type_id=1"
	. " AND a.value=c.ec"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1"
	. " GROUP by a.feature_id, cofactor";
    my $ecsth = $dbh->prepare($ecq);
    $ecsth->execute();
    while (my @rv = $ecsth->fetchrow_array) {
	if (defined $r->{$rv[0]}) {
	    push @{$r->{$rv[0]}->{'cofactor'}->{'ec'}}, $rv[1];
	}
    }
    
    my $xq = "select a.feature_id, source, prefix, accession"
	. " FROM feature_accessions a, seq_set_link l, seq_feat_mappings m, sequences s"
	. " WHERE set_id=$set_id"
	. " AND m.seq_id=l.seq_id"
	. " AND s.seq_id=m.seq_id AND s.iscurrent=1"
	. " AND a.feature_id=m.feature_id";
    my $xsth = $dbh->prepare($xq);
    $xsth->execute();
    while (my @rv = $xsth->fetchrow_array) {
	if (defined $r->{$rv[0]}) {
	    my $acc = $rv[2] ? "$rv[2]|$rv[3]" : "gnl|$rv[3]";
	    $r->{$rv[0]}->{'accessions'}->{$rv[1]} = $acc;
#	$r->{$acc} = \{$r->{$rv[0]} unless (defined $r->{$acc});
	}
    }
    return $r;
}

sub get_sequences_by_set_id {
    my $dbh = shift;
    my $set_id = shift;

    # get sequence
    my $sq = "select s.seq_id, sequence, seq_length, gc, category"
	. " FROM sequences s, seq_set_link l"
	. " WHERE set_id=$set_id"
	. " AND s.seq_id=l.seq_id"
	. " AND s.iscurrent=1"
	;
    my $seq_ref = $dbh->selectall_hashref($sq, 'seq_id');

    # get sequence annotations
    my $nq = "select a.seq_id, data_type, value"
	. " FROM sequence_annotations a, seq_set_link l"
	. " WHERE set_id=$set_id"
	. " AND a.seq_id=l.seq_id";
    my $nsth = $dbh->prepare($nq);
    $nsth->execute();
    while (my @rv = $nsth->fetchrow_array) {
	push @{$seq_ref->{$rv[0]}->{'annotation'}->{$rv[1]}}, $rv[2] if (defined ($seq_ref->{$rv[0]}));
    }

    # get sequence accessions
    my $xq = "select a.seq_id, seq_acc_source, seq_accession"
	. " FROM sequence_accessions a, seq_set_link l"
	. " WHERE set_id=$set_id"
	. " AND a.seq_id=l.seq_id";
    my $xsth = $dbh->prepare($xq);
    $xsth->execute();
    while (my @rv = $xsth->fetchrow_array) {
#	my $acc = $rv[1] ? "$rv[1]|$rv[2]" : $rv[2];
	my $acc = $rv[2];
	push @{$seq_ref->{$rv[0]}->{'accessions'}}, $acc if (defined ($seq_ref->{$rv[0]}));
    }

    return $seq_ref;
}

sub seq_id_to_set_name {
    my $dbh = shift;
    my $seq_id = shift;
    my $q = "SELECT s.name"
	. " FROM sequence_sets s, seq_set_link l"
	. " WHERE l.seq_id=$seq_id"
	. " AND s.set_id=l.set_id"
	. " ORDER BY s.set_id DESC";
    my $r = $dbh->selectall_arrayref($q);
    return $r->[0]->[0];
}

sub set_name_to_id {
    my $dbh =  shift;
    my $set_name = shift;
    my $q = "select set_id from sequence_sets where name=\"$set_name\"";
    my $r = $dbh->selectall_arrayref($q);
    if (@$r > 1) {
	warn "More than one id for this set name: $set_name\n";
	foreach my $row (@$r) {
	    print "\t$row->[0]";
	}
	print "\n";
	die();
    } elsif (@$r < 1) {
	warn "No set_id found for $set_name\n";
	return();
    } else {
	return $r->[0]->[0];
    }
}

sub set_id_to_set_name {
    my $dbh =  shift;
    my $set_id = shift;
    my $q = "select name from sequence_sets where set_id=$set_id";
    my $r = $dbh->selectall_arrayref($q);
    if (@$r > 1) {
	warn "More than one name for this set id: $set_id\n";
	foreach my $row (@$r) {
	    print "\t$row->[0]";
	}
	print "\n";
	die();
    } elsif (@$r < 1) {
	warn "No set_name found for $set_id\n";
	return();
    } else {
	return $r->[0]->[0];
    }
}

sub get_evidence_for_feature {
    my $dbh = shift;
    my $feature_id = shift;
    my $ev_type = shift;

    my $hmm_q = "select * from feature_evidence WHERE feature_id=$feature_id";
    if ($ev_type) { $hmm_q .= " AND ev_type=\"$ev_type\"" }

    # [0] feature_id
    # [1] feat_min
    # [2] feat_max
    # [3] program
    # [4] ev_type
    # [5] ev_accession
    # [6] ev_min
    # [7] ev_max
    # [8] ev_length
    # [9] score
    # [10] ev_date
    return $dbh->selectall_arrayref($hmm_q);
}

sub get_mainroles {
    my $dbh = shift;
    my $mrq = "select * from egad.mainrole";
    return $dbh->selectall_hashref($mrq, 'id');
}

sub get_subroles {
    my $dbh = shift;
    my $mrq = "select * from egad.subrole";
    return $dbh->selectall_hashref($mrq, 'id');
}

sub set_id_to_seq_ids {
    my $dbh = shift;
    my $set_id = shift;

    my $sq = "select s.seq_id"
	. " FROM sequences s, seq_set_link l"
	. " WHERE set_id=$set_id"
	. " AND s.seq_id=l.seq_id"
	. " AND s.iscurrent=1"
	. " ORDER BY s.seq_length desc";
    my $seq_ref = $dbh->selectcol_arrayref($sq);
    return $seq_ref;
}

sub seq_id_to_SeqObj {
    my $dbh = shift;
    my $seq_id = shift;
    
    my $sq = "select sequence"
	. " FROM sequences s"
	. " WHERE seq_id=$seq_id"
	;
    my @seqrow = $dbh->selectrow_array($sq);

    # get sequence accessions
    my $xq = "select seq_acc_source, seq_accession"
	. " FROM sequence_accessions a"
	. " WHERE seq_id=$seq_id"
	;
    my $xsth = $dbh->prepare($xq);
    $xsth->execute();
    my @rv = $xsth->fetchrow_array;

    my $SeqObj = Bio::Seq::RichSeq->new(-seq => $seqrow[0],
					-molecule => 'DNA',
					-display_name => $rv[1],
					-accession_number => $rv[1],
					-id => $rv[1],
					-primary_id => $rv[1],
					-alphabet => 'dna',
	);

    # get sequence annotations
    my $nq = "select data_type, value"
	. " FROM sequence_annotations a"
	. " WHERE seq_id=$seq_id"
	;
    my $nsth = $dbh->prepare($nq);
    $nsth->execute();
    while (my @rv = $nsth->fetchrow_array) {
	my $AnnObj;
	if ($rv[0] eq 'comment') {
	    $AnnObj = Bio::Annotation::Comment->new(-text => $rv[1]);
	} elsif ($rv[0] eq 'dblink') {
	    my ($db, $id) = split/\:/, $rv[1];
	    $AnnObj = Bio::Annotation::DBLink->new(-database => $db,
						   -primary_id => $id);
	} elsif ($rv[0] eq 'ontology_term') {
	    my ($ontology, $term) = split/\:/, $rv[1];
	    $AnnObj = Bio::Annotation::OntologyTerm->new(-term => $term,
							 -tagname => $ontology);
	} elsif ($rv[0] eq 'reference') {
	    my ($medline, $pubmed, $title, $authors, $location, $start, $end, $rp_line, $rg_line);
	    $AnnObj = Bio::Annotation::Reference->new(-medline => $medline,
						      -pubmed => $pubmed,
						      -title => $title,
						      -authors => $authors,
						      -location => $location,
						      -start => $start,
						      -end => $end,
						      -rp => $rp_line,
						      -rg => $rg_line);
	} elsif ($rv[0] eq 'description') {
	    $SeqObj->description($rv[1]);
	} else {
	    $AnnObj = Bio::Annotation::SimpleValue->new(-value => $rv[1],
							-tagname => $rv[0]);
	}
	$SeqObj->add_Annotation($AnnObj) if (defined $AnnObj);
    }
    if (! $SeqObj->description) {
	my $set_name = &seq_id_to_set_name($dbh, $seq_id);
	$SeqObj->description("$set_name scaffold " . $SeqObj->display_name);
    }
    return $SeqObj;
}

sub add_features_to_SeqObj {
    my $dbh = shift;
    my $seq_id = shift;
    my $SeqObj = shift;

    my $srcObj = new Bio::SeqFeature::Generic;
    $srcObj->primary_tag('source');
    $srcObj->start(1);
    $srcObj->end($SeqObj->length);
    $SeqObj->add_SeqFeature($srcObj);

    my $feat_id_a = &get_feature_ids_by_seq_id($dbh, $seq_id);
    my $feat_ref = &get_features_by_feature_id($dbh, @$feat_id_a);
    
    foreach my $feat_id (sort {$feat_ref->{$a}->{'location'}->{$seq_id}->{'feat_min'} <=>
				   $feat_ref->{$b}->{'location'}->{$seq_id}->{'feat_min'} }
			 keys %$feat_ref) {
	my $featr = $feat_ref->{$feat_id};

	my $FeatObj = new Bio::SeqFeature::Generic;
	$FeatObj->primary_tag($featr->{'feat_type'});
	$FeatObj->add_tag_value('translation', $featr->{'product'}) if ($featr->{'product'});

	### Valid qualifiers for CDS
#	/allele="text"
#	/artificial_location="[artificial_location_value]"
#	/citation=[number]
#	/codon_start=<1 or 2 or 3>
#       /db_xref="<database>:<identifier>"
#       /EC_number="text"
#       /exception="[exception_value]"
#       /experiment="[CATEGORY:]text"
#       /function="text"
#       /gene="text"
#       /gene_synonym="text"
#	/inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
#	/locus_tag="text" (single token)
#	/map="text"
#	/note="text"
#	/number=unquoted text (single token)
#	/old_locus_tag="text" (single token)
#	/operon="text"
#	/product="text"
#	/protein_id="<identifier>"
#	/pseudo
#	/pseudogene="TYPE"
#	/ribosomal_slippage
#	/standard_name="text"
#	/translation="text"
#	/transl_except=(pos:<base_range>,aa:<amino_acid>)
#	/transl_table =<integer>
#	/trans_splicing

	# annotations
	my @ranks = sort {$a<=>$b} keys %{$featr->{'annotation'}};
	# to only pull annotation from best rank:
	foreach my $source (keys %{$featr->{'annotation'}->{$ranks[0]}}) {
	    foreach my $ tag (keys %{$featr->{'annotation'}->{$ranks[0]}->{$source}}) {
		$FeatObj->add_tag_value($tag, @{$featr->{'annotation'}->{$ranks[0]}->{$source}->{$tag}});
	    }
	}

	# location
	my $LocObj;
	my $tcoords = $featr->{'location'}->{$seq_id}->{'translation_coords'};
        #location if translation_coords
	if($tcoords){
	    my $locfact = Bio::Factory::FTLocationFactory->new();
	    $LocObj = $locfact->from_string($tcoords);
#	    $LocObj= Bio::LocationI->new(-start => $loc->start(),
#					 -end => $loc->end(),
#					 -strand => $featr->{'location'}->{$seq_id}->{'strand'});
	} else{
	    #location if no translation_coords
	    my $start_fuzz = $featr->{'location'}->{$seq_id}->{'min_partial'} ? "BEFORE" : "EXACT";
	    my $end_fuzz = $featr->{'location'}->{$seq_id}->{'max_partial'} ? "AFTER" : "EXACT";
	    $LocObj=Bio::Location::Fuzzy->new(-start => $featr->{'location'}->{$seq_id}->{'feat_min'},
					      -end => $featr->{'location'}->{$seq_id}->{'feat_max'},
					      -strand => $featr->{'location'}->{$seq_id}->{'strand'},
					      -start_fuzz => $start_fuzz,
					      -end_fuzz => $end_fuzz);
	}
	$FeatObj->location($LocObj);
	$FeatObj->frame($featr->{'location'}->{$seq_id}->{'phase'}) if ($featr->{'location'}->{$seq_id}->{'phase'} =~ /\d/);
	
	# accessions as protein_ids
	my $locus = defined $featr->{'accessions'}->{'locus_tag'} ? $featr->{'accessions'}->{'locus_tag'}->[0] : 
	    defined $featr->{'accessions'}->{'PNNL'} ? $featr->{'accessions'}->{'PNNL'}->[0] : "";
	$FeatObj->add_tag_value('locus_tag', $locus) if ($locus);

	if ($featr->{'feat_type'} eq 'CDS') {
	    &add_Gene_feature($SeqObj, $featr, $LocObj);
	}

	# Finally, add FeatObj to SeqObj
	$SeqObj->add_SeqFeature($FeatObj);
    }
}

sub add_Gene_feature {
    my $SeqObj = shift;
    my $featr = shift;
    my $LocObj = shift;

    ### Valid qualifiers
    my @quals = ('allele',
		 'citation',
		 'db_xref',
		 'experiment',
		 'function',
		 'gene',
		 'gene_synonym',
		 'inference',
		 'locus_tag',
		 'map',
		 'note',
		 'old_locus_tag',
		 'operon',
		 'product',
		 'pseudo',
		 'pseudogene',
		 'phenotype',
		 'standard_name',
		 'trans_splicing');

    my $FeatObj = new Bio::SeqFeature::Generic;
    $FeatObj->primary_tag('gene');

    # location
    my $loc;
    if ($LocObj->isa('Bio::Location::SplitLocationI')) {
	$loc = Bio::Location::Simple->new(-start=>$LocObj->start,
					  -end=>$LocObj->end,
					  -strand=>$LocObj->strand);
    } else { $loc = $LocObj; }
    $FeatObj->location($loc);

    # annotations
    foreach my $tag (keys %{$featr->{'annotation'}}) {
	if (grep /^$tag$/, @quals) { 
	    $FeatObj->add_tag_value($tag, @{$featr->{'annotation'}->{$tag}});
	}
    }
    # accessions as locus_tag for now
    $FeatObj->add_tag_value('locus_tag', $featr->{'accessions'}->{'locus_tag'}->[0]) if (defined $featr->{'accessions'}->{'Locus_tag'});
    $SeqObj->add_SeqFeature($FeatObj);
}

sub load_seq_feat_mappings {
    my $dbh = shift;
    my $feat_id = shift;
    my $seq_id = shift;
    my $feato = shift;

    $feato->start =~ /(\D?)(\d+)/ || carp "What's with the start?: " . $feato->start . "\n";
    my ($start_partial, $start) = ($1,$2);
    my $min_partial = $start_partial eq "<" ||
	$feato->location->min_start != $feato->location->max_start ? 1 : 0;
    
    $feato->end =~ /(\D?)(\d+)/ || carp "What's with the end?: " . $feato->end . "\n";
    my ($end_partial, $end) = ($1,$2);
    my $max_partial = $end_partial eq ">" ||
	$feato->location->min_end != $feato->location->max_end ? 1 : 0;

    my $translation_coords = $feato->location->to_FTstring;
    
    my $prot = "";
    my $phase = 0;
    if (defined $feato->frame) {
	$phase = $feato->frame;
    } else {
	my @codon_start = $feato->annotation->get_Annotations('codon_start');
	if (defined $codon_start[0]) {
	    $phase = $codon_start[0]->value - 1;
	}
    }
    my $strand = $feato->strand;

    my $seq_feat_mappings_i = sprintf "INSERT seq_feat_mappings"
	. " (seq_id, feature_id, feat_min, feat_max, strand,"
	. " phase, min_partial, max_partial, translation_coords)"
	. " VALUES ($seq_id, $feat_id, $start, $end, '$strand',"
	. " '$phase', $min_partial, $max_partial, '$translation_coords')";
#    print STDERR $seq_feat_mappings_i, "\n";
    my $row2 = $dbh->do($seq_feat_mappings_i);
    if (! defined $row2) {
	warn "\n!!!Couldn't execute '$seq_feat_mappings_i'\n";
	return;
    }
}

sub is_current {
    my $dbh = shift;
    my $fid = shift;
    my $is_current = shift;
    if ($is_current ne "") {
	my $q = "update sequence_features set is_current=$is_current where feature_id=$fid";
	my $r = $dbh->do($q);
    } else {
	my $q = "SELECT f.is_current, m.seq_id"
	    . " FROM sequence_features f, seq_feat_mappings m, seq_set_link l, sequence_sets s"
	    . " WHERE f.feature_id=$fid"
	    . " AND m.feature_id=f.feature_id"
	    . " AND l.seq_id=m.seq_id"
	    . " AND s.set_id=l.set_id"
	    . " AND s.is_current=1";
	my $r = $dbh->selectall_arrayref($q);
	if (@$r > 1) { warn "feature $fid is current on more than one sequence in the current sequence set?\n" }
	my $current;
	foreach my $row (@$r) {
	    my ($c, $s) = @$row;
	    $current = $current == 0 ? $c : $current;
	}
	return $current;
    }
}

sub update_product {
    my $dbh = shift;
    my $fid = shift;
    my $product = shift;
    my $q = "UPDATE sequence_features set product=\"$product\" WHERE feature_id=$fid";
    my $r = $dbh->do($q);
    if (! defined $r) {
	warn "Couldn't execute '$q'\n" . $dbh->errstr . "\n";
    }
}

sub insert_feature_accessions {
    my $dbh = shift;
    my $fid = shift;
    my $acc = shift;
    my $source = shift;
    my $prefix = shift;

    my $feat_acc_i = "INSERT feature_accessions"
	. " (feature_id, source, prefix, accession)"
	. " VALUES ($fid, \"$source\", \"$prefix\", \"$acc\")" ;
    my $r = $dbh->do($feat_acc_i);

}

sub update_feature_mapping {
    my $dbh = shift;
    my $sid = shift;
    my $fid = shift;
    my $map_ref = shift;
    if (! $sid) { warn "Empty sequence_id provided.\n"; return }
    if (! $fid) { warn "Empty feature_id provided.\n"; return }
    if (! $map_ref) { warn "No updates provided.\n"; return }

    my $map_q = "SELECT 1 from seq_feat_mappings where feature_id=$fid and seq_id=$sid";
    my @r = $dbh->selectrow_array($map_q);
    if ($r[0] != 1) { warn "feature $fid does not currently map to sequence $sid\n" }
    else {
	my @list;
	while (my($key,$val) = each %$map_ref) {
	    push @list, "$key=\"$val\"" if (defined $val);
	}
	my $map_u = "UPDATE seq_feat_mappings set " . join(",", @list) . " WHERE seq_id=$sid and feature_id=$fid";
	print "$map_u\n";
	my $r = $dbh->do($map_u);
	if (!$r) { warn $DBI::errorstr, "\n"; }
    }
    return;
}

sub set_name_to_set_id {
    my $dbh = shift;
    my $set_name = shift;

    if (! $set_name) { croak "No set_name provided."; }

    my $q = "SELECT set_id from sequence_sets where name=\"$set_name\"";
    my $r = $dbh->selectcol_arrayref($q);
    if (@$r == 0) { carp "No set_id found for '$set_name'.\n"; return; }
    elsif (@$r > 1) { carp "Multiple set_ids for '$set_name'. Returning " . $r->[0] . "\n"; }
    return $r->[0];
}

our %PSEUDOKEY = ( 1 => 'singleFS',
		   2 => 'singleIntStop',
		   3 => 'singleIndel',
		   4 => 'multipleLesions',
		   5 => 'unspecified' );
1;
