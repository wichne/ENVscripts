sub ConnectToDb {
    my($dbtype,$server,$username,$password,$db) = @_;
    my($dbproc);
    $dbproc = (DBI->connect("dbi:$dbtype:server=$server",$username, $password));
    if ( !defined $dbproc ) {
	die "Cannot connect to Sybase server: $DBI::errstr\n";
    }
    
    $dbproc->do("use $db");
    return($dbproc);
}

sub do_sql{
    my($dbproc,$query) = @_;
    my(@results,@row);
    
    my($statementHandle) = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
        die "Cannot prepare statement: $DBI::errstr\n";
    }
    $statementHandle->execute() || die $query, $statementHandle->errstr;
    #Fetch the rows back from the SELECT statement;
    if($statementHandle->{syb_more_results} ne "") {
        while ( @row = $statementHandle->fetchrow() ) {
            push(@results,join("\t",@row));
        }
    }
    # Release the statement handle resources;
    $statementHandle->finish;
    return @results;
}

sub RunMod{
    my($dbproc,$query) = @_;
    my($statementHandle);
    my(@results,@row);
    print "$query\n" if($DEBUG);
    $statementHandle = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
	die "Cannot prepare statement: $DBI::errstr\n";
    }
    $statementHandle->execute() || die $query, $statementHandle->errstr;
    $statementHandle->finish;
}

sub single_sql {
    my($dbproc,$query) = @_;
    my(@data, @results, $check);
    @results = &do_sql($dbproc, $query);
    $check = @results;
    if ($check >1) {
        return $results[0];
    } elsif ($check < 1) {
    } elsif ($check ==1) {
        return $results[0];
    }
}

sub QueryDBtoHashRef {
    my($dbproc,$query,$DEBUG) = @_;
    my($statementHandle);

    print "query: $query\n" if($DEBUG);
    
    $statementHandle = $dbproc->prepare($query);
    if ( !defined $statementHandle) {
        print "Cannot prepare statement: $DBI::errstr\n";
    }
    $statementHandle->execute() || print "failed query: $query\n";

    $ref = $statementHandle->fetchall_arrayref({});
    #release the statement handle resources
    $statementHandle->finish;
    return($ref);
}

sub numerically {
    $a <=> $b;
}

sub GetAuthorFromId {
    my($id) = shift;
    my $name;
    my($query) = "select fname, mname, lname, suffix, initials "
	. "from gb_authors "
	    . "where id=$id";
    print "$query<br>\n" if($DEBUG == 1);
    my($y) = &do_sql($$dbhashref{'dbproc'},$query);
    my(@x) = split(/\t/,$y);
    return(@x);
}

sub CreateSelectFromGbType {
    my($gb_type_ref,$type,$default_id) = @_;
    my $select = "<select name=\"$type\">";
    foreach my $id( keys %$gb_type_ref) {
	if($$gb_type_ref{$id}->{'type'} eq $type) {
	    if ($id == $default_id) {
		$select .= "<option value=\"$id\" selected=1>$$gb_type_ref{$id}->{'type_desc'}";
		} 
	    else {
		$select .= "<option value= \"$id\">$$gb_type_ref{$id}->{'type_desc'}";
	    }
	}
    }
    $select .= "</select>";
    return($select);
}

sub CreateSelectFromGBFeatType {
    my($gb_feat_type_ref) = @_;
    my $isgene_select;
    my $istranslated_select;

    my $feat_type_select = "<select name =\"feat_type_id\">\n";
    my $gb_feat_type_select = "<select name =\"gb_feat_type_id\">\n";
    foreach my $id (sort keys %$gb_feat_type_ref) {
	my $feat_type = $gb_feat_type_ref->{$id}->{'feat_type'};
	my $gb_feat_type = $gb_feat_type_ref->{$id}->{'gb_feat_type'};
	$feat_type_select .= "<option value= \"$id\">$feat_type";
	$gb_feat_type_select .= "<option value= \"$id\">$gb_feat_type";
    }
    $feat_type_select .= "</select>";
    $gb_feat_type_select .= "</select>";

    $isgene_select = &CreateSelectFromYesNo("isgene");
    $istranslated_select = &CreateSelectFromYesNo("istranslated");
    return($feat_type_select,$gb_feat_type_select,$isgene_select,$istranslated_select);
}

sub CreateSelectFromYesNo {
    my($type) = @_;
    my $select = "<select name=\"$type\">";
    $select .= "<option value= \"1\">YES";
    $select .= "<option value= \"0\">NO";
    $select .= "</select>";
    return($select);
}
sub FormatNullString {
    my($value) = @_;
    $value =~ s/\A\s+//;
    $value =~ s/\s+\Z//;
    my($retval) = "NULL";
    if($value ne "") {
	$retval = "'$value'";
    }
    return($retval);
}
sub FormatNullNumber {
    my($value) = @_;
    $value =~ s/\s+//;
    if($value eq "") {
	$value = "NULL";
    }
    return $value;
}

sub LoadFeatType {
    my($feat_type, $gb_feat_type, $isgene, $istranslated) = @_;
    my $query = "insert gb_feat_type(feat_type,gb_feat_type,isgene,istranslated) \n"
	. "values("
	. &FormatNullString($feat_type).", "
	. &FormatNullString($gb_feat_type).", "
	    . "$isgene,$istranslated)\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
    return("");
}
sub UpdateFeatType {
    my($feat_type, $gb_feat_type, $isgene, $istranslated,$update_id) = @_;
    my $query = "update gb_feat_type \n"
	. "set feat_type = ". &FormatNullString($feat_type) . ", "
	. "gb_feat_type = ". &FormatNullString($gb_feat_type) . ", "
	. "isgene = $isgene, "
	. "istranslated = $istranslated "
	. "where id = $update_id\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub LoadFeatTypeLink {
    my($gb_feat_type_id,$gb_release_id) = @_;
    my $query = "insert gb_feat_type_link(gb_feat_type_id,gb_release_id) \n"
    . "values($gb_feat_type_id, $gb_release_id)\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub DeleteFeatTypeLink {
    my($gb_release_id) = @_;
    my $query = "delete gb_feat_type_link"
	. " where gb_release_id=$gb_release_id"
	;
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub DeleteAuthorsLink {
    my($gb_citations_id) = @_;
    my $query_delete = "delete gb_annotation..gb_author_link
                        where gb_citations_id = $gb_citations_id\n";
    &RunMod($$dbhashref{'dbproc'},$query_delete);
}

sub LoadAuthor {
    my($fname,$mname,$lname,$suffix,$initials) = @_;
    my $insert_query = "insert gb_authors(lname,fname,mname,suffix,initials)\n"
	. "values('$lname','$fname'," 
	.&FormatNullString($mname).","
	.&FormatNullString($suffix).","
	.&FormatNullString($initials).")\n";
    print "$insert_query\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$insert_query) if($DEBUG != 1);
}

sub LoadAuthorLink {
    my($gb_author_id,$gb_citations_id,$author_order) = @_;
    my $query = "insert gb_annotation..gb_author_link(gb_authors_id,gb_citations_id,author_order) \n"
    . "values($gb_author_id, $gb_citations_id, $author_order)\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub InsertCitation {
    my ($ref) = shift;
    my($gb_citation_status,$title,$journal,$vol,$pages,$id) = @_;
    
    my $query = "insert gb_citations (title, journal, vol, pages, gb_citation_status, date)"
	. "values (". &FormatNullString($ref->{'title'}) . ", " 
	. &FormatNullString($ref->{'journal'}) . ", "
	. &FormatNullNumber($ref->{'vol'}) . ", "
	. &FormatNullString($ref->{'pages'}) . ", "
	. "$ref->{'gb_citation_status'}, "
	. "getdate())";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);

    my $q2 = "select \@\@identity";
    my ($id) = &do_sql($$dbhashref{'dbproc'}, $q2);

    return ($id);
}

sub UpdateCitation {
    my ($ref) = shift;
    my($gb_citation_status,$title,$journal,$vol,$pages,$id) = @_;
    
    my $query = "update gb_citations "
	. "set gb_citation_status = $ref->{'gb_citation_status'}, "
	. "title = "     . &FormatNullString($ref->{'title'}) . ", "
	. "journal = "     . &FormatNullString($ref->{'journal'}) . ", "
	. "vol = "  . &FormatNullNumber($ref->{'vol'}) . ", "
	. "pages = "    . &FormatNullString($ref->{'pages'}) . " "
	. "where id = $ref->{'cit_id'}\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub LinkCitation {
    my $ref = shift;
    my $link_q = "insert gb_citations_link (gb_release_id, gb_citations_id, gb_citation_type)"
	. " values ($ref->{'rel_id'}, $ref->{'cit_id'}, $ref->{'gb_citation_type'})";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$link_q);
}

sub UpdateCitationLink {
    my $ref = shift;
    my $link_q = "update gb_citations_link"
	. " set gb_citation_type=$ref->{'gb_citation_type'})"
	. " where gb_citations_id=$ref->{'cit_id'}"
	. " and gb_release_id=$ref->{'rel_id'}";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$link_q);
}

sub UnlinkCitation {
    my ($ref, $id) = @_;
    my $link_query = "delete gb_citations_link"
	. " where gb_citations_id = $id"
	. " and gb_release_id = $ref->{'rel_id'}";
    print "$link_query<br>\n" if ($DEBUG);
    &RunMod($$dbhashref{'dbproc'},$link_query);
}

sub DeleteCitation {
    my $ref = shift;
    my $link_q = "select gb_release_id"
	. " from gb_citations_link"
        . " where gb_citations_id=$ref->{cit_id}";
    my @links = &do_sql($$dbhashref{'dbproc'}, $link_q);

    if (@links) { print "Cannot delete citation - it is still linked to " . @links . " releases.<p>\n"; }
    else {
	my $cit_d = "delete gb_citations where id=$ref->{cit_id}";
    }
}

sub GetAuthorIdFromName {
    my($fname,$lname,$mname,$initials,$suffix) = @_;
    my($id) = 0;
    my $query = "select id from gb_authors
                     where lname = '$lname'
                     and fname = '$fname'\n";
    $query .= "and mname = '$mname'\n" if($mname ne "");
    $query .= "and suffix = '$suffix'\n" if($suffix ne "");
    $query .= "and initials = '$initials'\n" if($initials ne "");
    print "query: $query\n" if($DEBUG == 1);
    my @y = &do_sql($$dbhashref{'dbproc'},$query);
    if(scalar(@y) == 0) {
	my $insert_query = "insert gb_authors(lname,fname,mname,suffix,initials)\n"
	    . "values('$lname','$fname'," 
		.&FormatNullString($mname).","
		    .&FormatNullString($suffix).","
			.&FormatNullString($initials).")\n";
	print "$insert_query\n" if($DEBUG == 1);
	&RunMod($$dbhashref{'dbproc'},$insert_query) if($DEBUG != 1);
	print "query: $query\n" if($DEBUG == 1);
	my @x = &do_sql($$dbhashref{'dbproc'},$query);
	$id = $x[0];
    }else{
	$id = $y[0];
    }
    die "Can not find an id for $fname $mname $lname $suffix $initials\n" if($id eq "" || $id == 0);
    return($id);
}
sub UpdateAuthor {
    my($fname,$lname,$mname,$initials,$suffix,$id) = @_;
    my $query = "update gb_authors "
	      . "set fname = " . &FormatNullString($fname) . ", "
	      . "lname = "     . &FormatNullString($lname) . ", "
	      . "mname = "     . &FormatNullString($mname) . ", "
	      . "initials = "  . &FormatNullString($initials) . ", "
	      . "suffix = "    . &FormatNullString($suffix) . " "
	      . "where id = $id\n";
    print "$query<br>\n" if($DEBUG == 1);
    &RunMod($$dbhashref{'dbproc'},$query);
}

sub GetGBFeatTypeLink {
    my ($rel_id) = shift;
    my(%gb_feat_type_ids);
    my($query) = "select gb_feat_type_id "
	. "from gb_feat_type_link"
	. " where gb_release_id=$rel_id";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	$gb_feat_type_ids{$y[$i]} = 1;
    }
    return(\%gb_feat_type_ids);
}

sub GetGBAuthorLink {
    my(%gb_author_link);
    $gb_author_link{'count'}= 0;
    my($query) = "select gb_authors_id, gb_citations_id, author_order "
	. "from gb_author_link";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	$gb_author_link{$gb_author_link{'count'}}->{'gb_authors_id'} = $x[0];
	$gb_author_link{$gb_author_link{'count'}}->{'gb_citations_id'} = $x[1];
	$gb_author_link{$gb_author_link{'count'}}->{'author_order'} = $x[2];
	$gb_author_link{'count'}++;
    }
    return(\%gb_author_link);
}

sub GetAuthorLinkFromCitationId {
    my ($cit_id) = shift;
    my $hashref = {};
    my($query) = "select gb_authors_id, author_order "
	. "from gb_author_link"
	    . " where gb_citations_id = $cit_id";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    foreach my $y (@y) {
	my(@x) = split(/\t/,$y);
	$hashref->{$x[1]} = $x[0];
    }
    return($hashref);
}

sub GetGBCitationsLinkFromReleaseId {
    my ($rel_id) = shift;
    my $cit_ids;
    my $query = "select gb_citations_id, gb_citation_type "
	. "from gb_citations_link"
	. " where gb_release_id=$rel_id";
    print "$query<br>\n" if($DEBUG == 1);
    $cit_ids = &QueryDBtoHashRef($$dbhashref{'dbproc'}, $query);
    return($cit_ids);
}

sub GetCitationIdsFromReleaseId {
    my($gb_citations_link_ref,$release_id,$use_id_ref) = @_;
    for(my $j = 0; $j < $$gb_citations_link_ref{'count'}; $j++) {
	my $this_release_id = $$gb_citations_link_ref{$j}->{'gb_release_id'}; 
	if($release_id == $this_release_id) {
	    my($citation_id) = $$gb_citations_link_ref{$j}->{'gb_citations_id'};
	    push(@$use_id_ref,$citation_id);
	}
    }
}

sub GetGBFeatType {
    my(%gb_feat_type);
    my($query) = "select id, feat_type, gb_feat_type, isgene, istranslated "
	. "from gb_feat_type";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	$gb_feat_type{$x[0]}->{'feat_type'} = $x[1];
	$gb_feat_type{$x[0]}->{'gb_feat_type'} = $x[2];
	$gb_feat_type{$x[0]}->{'isgene'} = $x[3];
	$gb_feat_type{$x[0]}->{'istranslated'} = $x[4];
    }
    return(\%gb_feat_type);
}

sub GetGBType {
    my(%gb_type);
    my($query) = "select id,type,value,type_desc "
	. "from gb_type";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	$gb_type{$x[0]}->{'type'} = $x[1];
	$gb_type{$x[0]}->{'value'} = $x[2];
	$gb_type{$x[0]}->{'type_desc'} = $x[3];
    }
    return(\%gb_type);
}

sub GetFeatTypeFromReleaseId {
    my ($rel_id) = @_;
    my(%feat_type);
    my $query = "select t.id, t.feat_type "
	. "from common..feat_type t, gb_feat_type_link l "
	    . "where gb_release_id=$rel_id "
		. "and t.id=l.gb_feat_type_id";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    foreach my $row (@y) {
	my @x = split (/\t/,$row);
	$feat_type{$x[0]}->{'feat_type'} = $x[1];
	$feat_type{$x[0]}->{'gb_feat_type'} = $x[2];
	$feat_type{$x[0]}->{'isgene'} = $x[3];
	$feat_type{$x[0]}->{'istranslated'} = $x[4];
    }
    return \%feat_type;
}

sub GetRelease {
    my(%release);
    $release{'count'} = 0;
    my($query) = "select id,asmbl_id,gb_accession,date,gb_genome_type,gb_molecule_type,gb_trans_table,title,comment,original_db,contact_author_id,assignby "
	. "from gb_release";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	$release{$release{'count'}}->{'id'} = $x[0];
	$release{$release{'count'}}->{'asmbl_id'} = $x[1];
	$release{$release{'count'}}->{'gb_accession'} = $x[2];
	$release{$release{'count'}}->{'date'} = $x[3];
	$release{$release{'count'}}->{'gb_genome_type'} = $x[4];
	$release{$release{'count'}}->{'gb_molecule_type'} = $x[5];
	$release{$release{'count'}}->{'gb_trans_table'} = $x[6];
	$release{$release{'count'}}->{'title'} = $x[7];
	$release{$release{'count'}}->{'comment'} = $x[8];
	$release{$release{'count'}}->{'original_db'} = $x[9];
	$release{$release{'count'}}->{'contact_author_id'} = $x[10];
	$release{$release{'count'}}->{'assignby'} = $x[11];
	$release{'count'}++;
    }
    return(\%release);
}

sub GetTopologyFromGbRelease {
    my($release_ref,$gb_type_ref,$gb_release_id) = @_;
    my($topology) = "circular";
    for(my $i = 0; $i < $$release_ref{'count'};$i++) {
	if($$release_ref{$i}->{'id'} == $gb_release_id) {
	    my $gb_molecule_type = $$release_ref{$i}->{'gb_molecule_type'};
	    if($gb_molecule_type == 435) {
		$topology = "linear";
	    }elsif($gb_molecule_type == 437) {
		$topology = "circular";
	    }
	}
    }
    return $topology;
}

sub GetReleaseFromReleaseId {
    my ($rel_id) = shift;
    my(%release);
    my($query) = "select id,asmbl_id,gb_accession,date,gb_genome_type,gb_molecule_type,gb_trans_table,title,comment,original_db,contact_author_id,assignby "
	. "from gb_release"
	    . " where id = $rel_id";
    print "$query<br>\n" if($DEBUG == 1);
    my ($y) = &do_sql($$dbhashref{'dbproc'},$query);
    my (@x) = split(/\t/,$y);
    $release{$rel_id}->{'asmbl_id'} = $x[1];
    $release{$rel_id}->{'gb_accession'} = $x[2];
    $release{$rel_id}->{'date'} = $x[3];
    $release{$rel_id}->{'gb_genome_type'} = $x[4];
    $release{$rel_id}->{'gb_molecule_type'} = $x[5];
    $release{$rel_id}->{'gb_trans_table'} = $x[6];
    $release{$rel_id}->{'title'} = $x[7];
    $release{$rel_id}->{'comment'} = $x[8];
    $release{$rel_id}->{'original_db'} = $x[9];
    $release{$rel_id}->{'contact_author_id'} = $x[10];
    $release{$rel_id}->{'assignby'} = $x[11];
    
    return(\%release);
}

sub GetCitationFromReleaseId {
    my $rel_id = shift;
    my(%citations);
    my($query) = "select c.id, c.date, c.gb_citation_status, c.title,"
	. " c.journal, c.vol, c.pages, l.gb_citation_type"
	. " from gb_citations c, gb_citations_link l"
	. " where l.gb_release_id = $rel_id"
	. " and l.gb_citations_id = c.id";
    print "$query<br>\n" if($DEBUG == 1);
    my $citations = &QueryDBtoHashRef($$dbhashref{'dbproc'},$query);
    return($citations);
}

sub GetCitationFromCitationId {
    my $id = shift;
    my(%citations);
    my($query) = "select c.date, c.gb_citation_status, c.title , c.journal, c.vol, c.pages"
	. " from gb_citations c"
	. " where c.id = $id";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	$citations{'date'} = $x[0];
	$citations{'gb_citation_status'} = $x[1];
	$citations{'title'} = $x[2];
	$citations{'journal'} = $x[3];
	$citations{'vol'} = $x[4];
	$citations{'pages'} = $x[5];
    }
    return(\%citations);
}

sub GetCitations {
    my(%citations);
    my($query) = "select id, date, gb_citation_status, title, journal, vol, pages "
	. "from gb_citations";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    foreach my $row (@y) {
	my(@x) = split(/\t/,$row);
	$citations{$x[0]}->{'date'} = $x[1];
	$citations{$x[0]}->{'gb_citation_status'} = $x[2];
	$citations{$x[0]}->{'title'} = $x[3];
	$citations{$x[0]}->{'journal'} = $x[4];
	$citations{$x[0]}->{'vol'} = $x[5];
	$citations{$x[0]}->{'pages'} = $x[6];
    }
    return(\%citations);
}

sub GetAuthorHashFromReleaseId {
    my($gb_citations_link_ref,$gb_author_link_ref,$gb_authors_ref,$gb_release_id_i) = @_;
    my(%authors_hash);
    $authors_hash{'count'} = 0;
    for(my $i = 0; $i < $$gb_citations_link_ref{'count'};$i++) {
	next if($$gb_citations_link_ref{$i}->{'gb_release_id'} != $gb_release_id_i);
	my $citation_id = $$gb_citations_link_ref{$i}->{'gb_citations_id'};
	foreach my $j (sort_multihashkeys_by_value($gb_author_link_ref,"author_order")) {
	    next if($j eq "count");
	    next if($$gb_author_link_ref{$j}->{'gb_citations_id'} != $citation_id);
	    my $author_order = $$gb_author_link_ref{$j}->{'author_order'};
	    my $gb_authors_id = $$gb_author_link_ref{$j}->{'gb_authors_id'};
	    for(my $k = 0; $k < $$gb_authors_ref{'count'};$k++) {
		next if($$gb_authors_ref{$k}->{'id'} != $gb_authors_id);
		my $id          = $$gb_authors_ref{$k}->{'id'};
		my $lname       = $$gb_authors_ref{$k}->{'lname'};
		my $fname       = $$gb_authors_ref{$k}->{'fname'};
		my $mname       = $$gb_authors_ref{$k}->{'mname'};
		my $suffix      = $$gb_authors_ref{$k}->{'suffix'};
		my $initials    = $$gb_authors_ref{$k}->{'inititals'};
		my $author_name = "$fname $mname $lname $suffix $initials";
		$authors_hash{$authors_hash{'count'}}->{'author_name'} = $author_name;
		$authors_hash{$authors_hash{'count'}}->{'author_order'} = $author_order;
		$authors_hash{$authors_hash{'count'}}->{'gb_authors_id'} = $id;
		$authors_hash{$authors_hash{'count'}}->{'gb_release_id'} = $gb_release_id_i;
		$authors_hash{$authors_hash{'count'}}->{'lname'} = $lname;
		$authors_hash{$authors_hash{'count'}}->{'fname'} = $fname;
		$authors_hash{$authors_hash{'count'}}->{'mname'} = $mname;
		$authors_hash{$authors_hash{'count'}}->{'suffix'} = $suffix;
		$authors_hash{$authors_hash{'count'}}->{'inititals'} = $inititals;
		$authors_hash{'count'}++;
	    }
	}
	last;
    }
    return(\%authors_hash);
}

sub GetAuthorsFromCitationId {
    my (@cit_ids) = @_;
    my %authors;
    foreach my $id(@cit_ids) {
	my $query = "select a.id, fname, mname, lname, suffix, initials"
	    . " from gb_authors a, gb_author_link l"
		. " where l.gb_citations_id=$id"
		    . " and a.id=l.gb_authors_id";
	print "$query<br>\n" if($DEBUG == 1);
	my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
	foreach my $y(@y) {
	    my @x = split /\t/, $y;
	    if (!defined $authors{$x[0]}) {
		$authors{$x[0]}->{'fname'} = $x[1];
		$authors{$x[0]}->{'mname'} = $x[2];
		$authors{$x[0]}->{'lname'} = $x[3];
		$authors{$x[0]}->{'suffix'} = $x[4];
		$authors{$x[0]}->{'initials'} = $x[5];
	    }
	}
    }
    return (\%authors);
}

sub SplitAuthorName {
    my($author_name) = @_;
    my($fname,$mname,$lname,$suffix,$initials) = split(/ /,$author_name);
    return($lname,$fname,$mname,$suffix,$initials);
}

sub GetAuthors {
    my(%authors);
    my($query) = "select id, lname, fname, mname, suffix, initials "
	. "from gb_authors "
	    . "order by lname";
    print "$query<br>\n" if($DEBUG == 1);
    my(@y) = &do_sql($$dbhashref{'dbproc'},$query);
    for(my $i = 0; $i < scalar(@y); $i++) {
	my(@x) = split(/\t/,$y[$i]);
	next if($x[1] eq "" && $x[2] eq "");
	$authors{$x[0]}->{'lname'} = $x[1];
	$authors{$x[0]}->{'fname'} = $x[2];
	$authors{$x[0]}->{'mname'} = $x[3];
	$authors{$x[0]}->{'suffix'} = $x[4];
	$authors{$x[0]}->{'inititals'} = $x[5];
    }
    return(\%authors);
}

sub GetAuthorsString {
    my($cit_ref, $cit_id, $authors_ref) = @_;
    $auth_ord_ref = $cit_ref->{$cit_id}->{'authors'};
    my($authors_string1) = "        std {\n"; #13 spaces
    my($authors_string2) = "          std {\n"; #24 spaces
    my($authors_string3) = "          std {\n"; #24 spaces
    my($space1) = "                "; #21 spaces
    my($space2) = "                 "; #32 spaces
    foreach my $i ( sort numerically keys %$auth_ord_ref) {
	my $first = $authors_ref->{$auth_ord_ref->{$i}}->{'fname'};
	my $middle = $authors_ref->{$auth_ord_ref->{$i}}->{'mname'};
	my $last = $authors_ref->{$auth_ord_ref->{$i}}->{'lname'};
	my $suffix = $authors_ref->{$auth_ord_ref->{$i}}->{'suffix'};
	my $initials = $authors_ref->{$auth_ord_ref->{$i}}->{'initials'};
	
	my (@authors) = ();
	push(@authors,"last \"$last\"") if($last);
	push(@authors,"first \"$first\"") if($first);
	push(@authors,"middle \"$middle\"") if($middle);
	push(@authors,"inititals \"$inititals\"") if($initials);
	push(@authors,"suffix \"$suffix\"") if($suffix);
	$authors_string1 .= "          {\n";  #14 spaces
	$authors_string1 .= "            name\n"; #17 spaces
	$authors_string1 .= "              name {\n"; #19 spaces
	$authors_string2 .= "            {\n";  #25 spaces
	$authors_string2 .= "              name\n"; #28 spaces
	$authors_string2 .= "                name {\n"; #30 spaces
	$authors_string3 .= "            {\n";  #25 spaces
	$authors_string3 .= "              name\n"; #28 spaces
	$authors_string3 .= "                name {\n"; #30 spaces
	for(my $j = 0; $j < scalar(@authors);$j++) {
	    $authors_string1 .= $space1 . $authors[$j] . " ,\n" if($j < scalar(@authors) - 1);
	    $authors_string1 .= $space1 . $authors[$j] . " } } ,\n" if($j == scalar(@authors) - 1);
	    $authors_string2 .= $space2 . $authors[$j] . " ,\n" if($j < scalar(@authors) - 1);
	    $authors_string2 .= $space2 . $authors[$j] . " } } ,\n" if($j == scalar(@authors) - 1);
	    $authors_string3 .= $space2 . $authors[$j] . " ,\n" if($j < scalar(@authors) - 1);
	    $authors_string3 .= $space2 . $authors[$j] . " } } ,\n" if($j == scalar(@authors) - 1);
	}
    }
    $authors_string1 =~ s/ \,\n\Z/ \} \,/;
    $authors_string2 =~ s/ \,\n\Z/ \} \} \,/;
    $authors_string3 =~ s/ \,\n\Z/ \} \,/;
    return($authors_string1,$authors_string2,$authors_string3);
}

sub GetLineageString {
    my($lineage_ref) = @_;
    my($lineage);
    for(my $i = 0; $i < $$lineage_ref{'count'};$i++) {
	$lineage .= "$$lineage_ref{$i}; ";
    }
    $lineage =~ s/\; $//;
    return $lineage;
}

sub GetTemplateFile {
    my($file,$contact_last,$contact_first,$citations_ref,$authors_ref,$date_year,$date_month,$date_day,$lit_status,$acc_status,$seq_len,$topology,$accession,$genetic_code,$lineage_ref,$lit_comment) = @_;
    my $org_name = $lineage_ref->{'organism_name'};
    
    #Which citation info is supposed to be used for this? The 'citation for accession reference' or Publication Descriptor
    my ($cit_id, $title, $journal, $vol, $pages) = &GetPublicationDescriptor($citations_ref);
    my ($authors_string1, $authors_string2, $authors_string3) = &GetAuthorsString($citations_ref, $cit_id, $authors_ref);
    my $lineage = &GetLineageString($lineage_ref);
    my($accession2) = $accession;
    $accession2 = "NULL" if($accession eq "");
    my($tigr_contact) = "owhite\@tigr.org";
    #Seq-submit
    #Submit-block
    #Put this back in when I can find out the number codes
    #subtype $acc_status ,
    #            prepub $lit_status } } } } }
    
    my $template_file = "Submit-block ::= {
  contact {
    contact {
      name
        name {
          last \"$contact_last\" ,
          first \"$contact_first\" } ,
      affil
        std {
          affil \"The Institute for Genomic Research\" ,
          city \"Rockville\" ,
          sub \"MD\" ,
          country \"USA\" ,
          street \"9712 Medical Center Dr\" ,
          email \"$tigr_contact\" ,
          fax \"301-828-0208\" ,
          phone \"301-838-0200\" ,
          postal-code \"20850\" } } } ,
  cit {
    authors {
      names
$authors_string1
      affil
        std {
          affil \"The Institute for Genomic Research\" ,
          city \"Rockville\" ,
          sub \"MD\" ,
          country \"USA\" ,
          street \"9712 Medical Center Dr\" ,
          postal-code \"20850\" } } ,
    medium email ,
    date
      std {
        year $date_year ,
        month $date_month ,
        day $date_day } ,
    descr \"$acc_status\" } ,
  reldate
    std {
      year $date_year ,
      month $date_month ,
      day $date_day } }
Seqdesc ::= pub {
  pub {
    article {
      title {
        name \"$title\" } ,
      authors {
        names
$authors_string2
      from
        journal {
          title {
            jta \"$journal\" } ,
          imp {
            date
              std {
                year $date_year ,
                month $date_month ,
                day $date_day } ,
            volume \"$vol\" ,
            pages \"$pages\" ,
            prepub submitted } } } } }
Seq-submit ::= {
  sub {
    contact {
      contact {
        name
          name {
            last \"$contact_last\" ,
            first \"$contact_first\" } ,
        affil
          std {
            affil \"The Institute for Genomic Research\" ,
            city \"Rockville\" ,
            sub \"MD\" ,
            country \"USA\" ,
            street \"9712 Medical Center Dr\" ,
            email \"$tigr_contact\" ,
            fax \"301-828-0208\" ,
            phone \"301-838-0200\" ,
            postal-code \"20850\" } } } ,
    cit {
      authors {
        names
$authors_string1
        affil
          std {
            affil \"The Institute for Genomic Research\" ,
            div \"general\" ,
            city \"Rockville\" ,
            sub \"MD\" ,
            country \"USA\" ,
            street \"9712 Medical Center Dr\" ,
            postal-code \"20850\" } } ,
      date
        std {
          year $date_year ,
          month $date_month ,
          day $date_day } } ,
    subtype 1 ,
    tool \"Sequin 3.30\" } ,
  data
    entrys {
      seq {
        id {
          local
            str \"$accession2\" } ,
        descr {
          title \"$title\" ,
          molinfo {
            biomol genomic } ,
          create-date
            std {
              year $date_year ,
              month $date_month ,
              day $date_day } ,
          source {
            genome genomic ,
            org {
              taxname \"$org_name\" ,
              orgname {
                lineage \"$lineage\" ,
                gcode $genetic_code ,
                div \"BCT\" } } } } ,
        inst {
          repr raw ,
          mol dna ,
          length $seq_len ,
          seq-data
            ncbi2na ''H } } } }\n";

 open(OUT,">$file") || die "Can not open file $file for writing\n";
 print OUT $template_file;
 close(OUT);
}

sub GetPublicationDescriptor {
    my ($cit_ref) = @_;
    foreach my $id (sort numerically keys %$cit_ref) {
	if ($cit_ref->{$id}->{'gb_citation_type'} != 433 ) {
	    return ($id, $cit_ref->{$id}->{'title'}, $cit_ref->{$id}->{'journal'},$cit_ref->{$id}->{'vol'}, $cit_ref->{$id}->{'pages'});
	}
    } 
}

sub CreateEditAuthorButton {
    my $button = "<input type=button name=edit_gb_author value=\"Add/Edit Author\" onclick=\"window.open('add_gb_authors.cgi?user=$$dbhashref{'user'}&password=$$dbhashref{'password'}')\">";
    return($button);
}

sub CreateDoneButton { 
    $html = "<p><input type=button name=done value=\"Done\" onclick=\"window.opener.location.reload();window.close();\">";
    return($html);
}
1;


