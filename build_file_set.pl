#!/usr/bin/perl

# This script is used to generate flatfiles of either
# all sequence_sets in the ENV database chosen,
# or a specified list
use lib $ENV{SCRIPTS};
use ENV;
use Getopt::Long;

print "USAGE: build_file_set.pl -name @set_names -id @set_ids -DB database -path dest_path

If no set list is provided by name or id, all current sets will be dumped.
If no path is specified, files written to current directory.
";
&GetOptions("name=s@" => \@set_names,
	    "id=i@" => \@set_ids,
	    "DB=s" => \$db,
	    "path=s" => \$dest_path);

my $dbh = &connect({ 'u' => 'access',
		     'p' => 'access',
		     'D' => $db } );

@set_names = split(/,/,join(',',@set_names));
@set_ids = split(/,/,join(',',@set_ids));

# We're going to build the following files:
#   contig/scaffold fasta
#   peptide sequence fasta
#   gene sequence fasta
#   GFF3
#   tab-delimited annotation file

foreach my $id(@set_ids) {
    my $name = set_id_to_set_name($dbh, $id);
    if ($name) { push @set_names, $name; }
}

if (! @set_names) {
    print "Dumping info for all sequence_sets in the database (Hope that's what you wanted to do...)\n";
    $set_ref = get_set_names($dbh);
    while (my($id, $ref) = each %$set_ref) {
	push @set_ids, $id;
	push @set_names, $ref->{'name'};
    }
}
foreach my $name (@set_names) {
    my $path = $dest_path ? $dest_path : ".";
    if (! -e $path) {
	mkdir $path, 0777 || die "Can't make output path directory $path\n";
    }
    if (! -w $path) {
	die "Can't write to specified output path $path\n";
    }

    
    $contig_cmd = "$ENV{SCRIPTS}/set_to_fasta.pl -D $db -n $name -o $path/$name.fasta";
    $pep_cmd    = "$ENV{SCRIPTS}/set_to_pep.pl -D $db -n $name -o $path/$name.pep";
    $seq_cmd    = "$ENV{SCRIPTS}/set_to_seq.pl -D $db -n $name -o $path/$name.seq -A locus_tag -a";
    $gff_cmd    = "$ENV{SCRIPTS}/set_to_gff.pl -D $db -n $name -o $path/$name.gff";
    $gbf_cmd    = "$ENV{SCRIPTS}/set_to_gbf.pl -D $db -n $name -o $path/$name.gbf";
    $csv_cmd    = "$ENV{SCRIPTS}/set_to_annotation_table.pl -D $db -n $name -o $path/$name.csv";

    print "$contig_cmd\n";
    my $r = system($contig_cmd);
    if ($r) { warn "WARNING: command '$contig_cmd' returned with $r: $!\n"; }

    print "$pep_cmd\n";
    $r = system($pep_cmd);
    if ($r) { warn "WARNING: command '$pep_cmd' returned with $r: $!\n"; }

    print "$seq_cmd\n";
    $r = system($seq_cmd);
    if ($r) { warn "WARNING: command '$seq_cmd' returned with $r: $!\n"; }

    print "$gff_cmd\n";
    $r = system($gff_cmd);
    if ($r) { warn "WARNING: command '$gff_cmd' returned with $r: $!\n"; }

    print "$gbf_cmd\n";
    $r = system($gbf_cmd);
    if ($r) { warn "WARNING: command '$gbf_cmd' returned with $r: $!\n"; }

    print "$csv_cmd\n";
    $r = system($csv_cmd);
    if ($r) { warn "WARNING: command '$csv_cmd' returned with $r: $!\n"; }


}

exit();
