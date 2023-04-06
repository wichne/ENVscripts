#!/usr/bin/perl
use DBI;
use Getopt::Std;
use strict;
#use lib $ENV{SCRIPTS};
#use ENV;
use Bio::FeatureIO;

# to use this program
# need to provide -c mcl_clusterfile? < *.ortho (these files are the output from inparanoid?)
my %arg;
&getopts('c:',\%arg);

my %GENECLUST;
my $clustercount=0;
open my $clusterIN, $arg{'c'} if ($arg{'c'});
# read in the mcl file
while (my $line = <$clusterIN>) {
    chomp $line;
    my @m = split/\t/, $line;
    $clustercount++;
    # assign current index to each accession
    $GENECLUST{$_} = $clustercount for @m;
}

my $GENES = {};
my %LOOKUP;
my %SYN;
my %SYNGROUP;
my $syngroupcount;
my %GROUP;
my %INFO;

# read in the bbh filenames
while (my $file = <STDIN>) {
    chomp $file;
    my %BBH;
    my ($set1, $set2);
    # parse the dataset names
    if ($file =~ /([^\.]+)\.([^\.]+)\.ortho/) {
	$set1 = $1;
	$set2 = $2;
    } else { die "Set names can't be determined from filename $file\n";}

    # read the bbhs into a structure
    open (my $bbh_in, $file) or die "Can't open $file: $!\n";

    # !! for the multiparanoid format of bbh reporting:
    # lines are in pairs, so read in the first line
    while (my $line1 = <$bbh_in>) {
	my @f1 = split/\s+/, $line1;
	my @acc1 = split/\|/, $f1[4];
	# now read in the second line
	my $line2 = <$bbh_in>;
	my @f2 = split/\s+/, $line2;
	if ($f1[0] != $f2[0]) { die "Problem with line pairing in bbh file $file:\n$line1$line2"; }
	my @acc2 = split/\|/, $f2[4];

    # !! - okay, if the bbh are in the 'basic' report format use the following code:
#    while (my $line = <$bbh_in>) {
#	chomp $line;
#	my @acc = split(/\s+/, $line);
#	if (@acc > 2) { die "Which bbh file format is this now? I'm expecting 'acc1\tacc2'. Instead I get: '$line'\n";}
#	my @acc1 = split(/\|/, $acc[0]);
#	my @acc2 = split(/\|/, $acc[1]);

	# assign the bbh
	my $acc1 = @acc1 > 1 ? $acc1[1] : $acc1[0]; # pick an accession
	my $acc2 = @acc2 > 1 ? $acc2[1] : $acc2[0];
	if (defined $BBH{$acc1}) { die "How can $acc1 already have a bbh assigned ($BBH{$acc1}) when I get to $acc2?\n";}
	if (defined $BBH{$acc2}) { die "How can $acc2 already have a bbh assigned ($BBH{$acc2}) when I get to $acc1?\n";}
	$BBH{$acc1} = $acc2;
	$BBH{$acc2} = $acc1;
    }
    
    # get sequences and genes from the sets and make ordered structures
    foreach my $set ($set1, $set2) {
	if (-e "$set.gff" && ! defined $INFO{$set}) {
	    my %D;
	    my $in  = Bio::FeatureIO->new(-file => $gff_file , -format => 'GFF');
	    while (my $featureo = $in->next_feature()) {
		if ($featureo->primary_tag eq "CDS") {
		    my $protid = $featureo->primary_id();
		    if ($protid =~ /cds-/) { $protid =~ s/cds-// }
		    $D{$protid} = {'posn' => $featureo->start(),
				       'dir' => $featureo->strand(),
				       'seq_id' => $featureo->seq_id() }

		}
	    }
#	    open my $info, "$set.gff";
#	    while (my $line = <$info>) {
#		chomp $line;
#		my ($feat_acc, $lo, $dir, $seq_id) = split/\s+/, $line;
#		my @acc = split/\|/, $feat_acc;
#		my $acc = @acc > 1 ? $acc[1] : $acc[0];
#		$D{$acc} = { 'posn' => $lo,
#			     'dir' => $dir,
#			     'seq_id' => $seq_id};
#	    }

	    # %D just holds info for the current set. Sort by seq_acc then by position.
	    foreach my $acc (sort {$D{$a}->{seq_id} cmp $D{$b}->{seq_id} || $D{$a}->{posn} <=> $D{$b}->{posn}} keys %D) {
		# for each seq_id, create an array of feature_ids ordered by position 
		push @{$GENES->{$set}->{$D{$acc}->{seq_id}}}, $acc;
		$LOOKUP{$acc} = {'set' => $set,
				 'seq' => $D{$acc}->{seq_id},
				 'idx' => $#{$GENES->{$set}->{$D{$acc}->{seq_id}}}};
	    }
	}
    }

    
    while (my ($g, $h) = each %BBH) {
	# only do comparison in one direction
	unless (($g cmp $h) < 0) { next }
	# see if the upstream and downstream genes are bbhs with 
	# the bbh's upstream and downstream genes
	my $usg1 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}+1];
	my $dsg1 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}-1];
	my $ush1 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}+1];
	my $dsh1 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}-1];
	my $usg2 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}+2];
	my $dsg2 = $GENES->{$LOOKUP{$g}->{'set'}}->{$LOOKUP{$g}->{'seq'}}->[$LOOKUP{$g}->{'idx'}-2];
	my $ush2 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}+2];
	my $dsh2 = $GENES->{$LOOKUP{$h}->{'set'}}->{$LOOKUP{$h}->{'seq'}}->[$LOOKUP{$h}->{'idx'}-2];
#	print "$usg\t$g\t$dsg  <=>  $ush\t$h\t$dsh\n$BBH{$usg}\t$BBH{$g}\t$BBH{$dsg}  <=>  $BBH{$ush}\t$BBH{$h}\t$BBH{$dsh}\n\n" if (($BBH{$usg} && ($BBH{$usg} eq $ush || $BBH{$usg} eq $dsh)) || ($BBH{$dsg} && ($BBH{$dsg} eq $ush || $BBH{$dsg} eq $dsh)));
	my $syn_score = 0;
	foreach my $sg ($usg1, $dsg1, $usg2, $dsg2) {
	    foreach my $sh ($ush1, $dsh1, $ush2, $dsh2) {
		$syn_score += 2 if ($BBH{$sg} && ($BBH{$sg} eq $sh));
		$syn_score += 1 if ($GENECLUST{$sg} && ($GENECLUST{$sg} eq $GENECLUST{$sh}));
	    }
	}

	if ($syn_score > 1) {
	    if (defined $SYNGROUP{$g} && defined $SYNGROUP{$h}){
		if ($SYNGROUP{$g} == $SYNGROUP{$h}) {
		} else {
		    print STDERR "$g and $h are already in synteny groups $SYNGROUP{$g} and $SYNGROUP{$h}\n";
		    print STDERR "\tTrying to merge groups $SYNGROUP{$g} and $SYNGROUP{$h}\n";
		    my $deletegroup = $SYNGROUP{$h};
		    my $conflicts;
		    foreach my $set (keys %{$GROUP{$deletegroup}}) {
			if (defined $GROUP{$SYNGROUP{$g}}->{$set}) {
			    print STDERR "\tCONFLICT: group $SYNGROUP{$g} already has a member for $set: $GROUP{$SYNGROUP{$g}}->{$set}\n";
			    $conflicts++;
			}
		    }
		    if ($conflicts == 0) {
			foreach my $set (keys %{$GROUP{$deletegroup}}) {
			    $SYNGROUP{$GROUP{$deletegroup}->{$set}} = $SYNGROUP{$g};
			    $GROUP{$SYNGROUP{$g}}->{$set} = $GROUP{$deletegroup}->{$set};
			}
			delete $GROUP{$deletegroup};
		    } else { 
			print STDERR "\tMerge failed.\n";
		    }
		}
	    } elsif (defined $SYNGROUP{$g}) {
		$SYNGROUP{$h} = $SYNGROUP{$g};
		$GROUP{$SYNGROUP{$h}}->{$LOOKUP{$h}->{'set'}} = $h;
	    } elsif (defined $SYNGROUP{$h}) {
		$SYNGROUP{$g} = $SYNGROUP{$h};
		$GROUP{$SYNGROUP{$g}}->{$LOOKUP{$g}->{'set'}} = $g;
	    } else {
		$syngroupcount++;
		$SYNGROUP{$g} = $syngroupcount;
		$SYNGROUP{$h} = $syngroupcount;
		$GROUP{$SYNGROUP{$g}}->{$LOOKUP{$g}->{'set'}} = $g;
		$GROUP{$SYNGROUP{$h}}->{$LOOKUP{$h}->{'set'}} = $h;
	    }
	}
    }
}

# print a header row
print "Id\t" . join("\t", sort keys %$GENES) . "\n";

foreach my $i (sort {$a<=>$b} keys %GROUP) {
    print "$i";
    foreach my $g(sort keys %$GENES) {
	print "\t" . $GROUP{$i}->{$g};
    }
#    print "$g\t";
#    &explode_print(\%SYN, $g);
    print "\n";
}

sub explode_print {
    my $href = shift;
    my $g = shift;
    foreach my $h (@{$href->{$g}->{2}}) {
	print "$h\t";
	&explode_print($href, $h);
    }
#    print "(";
#    foreach my $h (@{$href->{$g}->{1}}) {
#	&explode_print($href, $h);
#    }
#    print ")\n";
#    print "\n";
}
