#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Std;

my %arg;
&getopts('D:u:p:', \%arg);
my $dbh=&connect(\%arg);

while (my $l=<STDIN>) {
    chomp $l;
    my ($fid, $sid1, $sid2) = split/\s+/, $l;
    my $fref = &get_features_by_feature_id($dbh, $fid);
    my $feat_type = $fref->{$fid}->{'feat_type'};

    # insert a new feature and get the feature_id
    my $feature_i = "INSERT sequence_features"
	. " (feat_type_id, feat_type, inserted_by, date_inserted)"
	. " SELECT id, \"$feat_type\", USER(), NOW()"
	. " FROM INSDC.feature_key"
	. " WHERE feature_key = \"$feat_type\"";

    my $ins_flag = $dbh->do($feature_i);
    my $new_fid = $dbh->last_insert_id("%", "%", "", "");

    # copy over the annotation
    my ($k) = sort {$a<=>$b} keys %{$fref->{$fid}->{'annotation'}};
    while (my ($src, $qref) = each %{$fref->{$fid}->{'annotation'}->{$k}}) {
	foreach my $qual (keys  %{$fref->{$fid}->{'annotation'}->{$k}->{$src}}) {
	    my $val = $fref->{$fid}->{'annotation'}->{$k}->{$src}->{$qual}->[0];
	   
	    my $feat_ann_i = "INSERT into feature_annotations"
		. " (feature_id, data_type_id, value, source)"
		. " SELECT $new_fid, d.id, \"$val\", \"$src\""
		. " FROM INSDC.qualifier d where d.qualifier = \"$qual\"";
	    print STDERR $feat_ann_i;
	    my $ann_flag = $dbh->do($feat_ann_i);
	}
    }

    # update the seq_feat_mappings row
    my $sfm_u = "UPDATE seq_feat_mappings set feature_id=$new_fid"
	. " WHERE feature_id=$fid and seq_id = $sid2";
    print STDERR $sfm_u;
    my $sfm_flag = $dbh->do($sfm_u);
}
