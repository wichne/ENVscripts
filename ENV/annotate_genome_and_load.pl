#!/usr/bin/perl -w

use Getopt::Long;
use lib $ENV{ENVSCRIPTS};
use ENV;
use strict;

# need to start with a fasta file of contigs/scaffolds that have stable contig identifiers
my $PROKKA = "prokka";
my $HMMSEARCH = "hmmsearch";
my $TIGRFAMDB = "/projects/db/hmms/TIGRFAM/TIGRFAMs_15.0_HMM.LIB";
my $PFAMDB = "/projects/db/hmms/Pfam/Pfam-A.hmm";
my $KOFAMSCAN = "/projects/bifx/kofam_scan-1.3.0/exec_annotation";
my $KOFAMDB = "/projects/db/hmms/KOfam/KOFAM.hmm";
my $DBCAN = "run_dbcan.py";
my $DBCANDB = "/projects/db/hmms/dbCAN/V9/";
my $TMHMM = "/projects/bifx/tmhmm-2.0c/bin/tmhmm";
my $PSORT = "psort";

my($genome_fasta, $protein_faa, $organism, $arch, $gram, $cpus, $locus_prefix,
   $file_prefix, $db, $user, $pswd, $set_id, $help);
&GetOptions("fasta=s" => \$genome_fasta,
	    "protein=s" => \$protein_faa,
	    "organism=s" => \$organism,
	    "archaea" => \$arch,
	    "gram=s" => \$gram,
	    "cpus=i" => \$cpus,
	    "locus_prefix=s" => \$locus_prefix,
	    "prefix=s" => \$file_prefix,
	    "db=s" => \$db,
	    "user=s" => \$user,
	    "pswd=s" => \$pswd,
	    "set_id=i" => \$set_id,
	    "help" => \$help,
#	    "prokka" => \$run_prokka
    );
if ($help) {
    print STDERR "
    annotate_genome_and_load.pl --fasta /path/to/genome_fasta --protein /path/to/protein_faa --organism organism --archaea --gram +|- --cpus n --prefix file_prefix --locus_prefix locus_prefix --db database --user username --pswd password --set-id n --help\n";
    exit;
}
if (! $user) { $user = $ENV{USER} }

if ($organism) {
    my ($genus, $species, $strain) = split(/\s+/, $organism, 3);
}

if ($gram) {
    if ($gram eq "-" || $gram eq "n") { $gram = "neg"}
    elsif ($gram eq "+" || $gram eq "p") { $gram = "pos" }
    if ($gram ne "pos" && $gram ne "neg") { die "Can't interpret gram value '$gram'. Please use 'pos' or 'neg'\n"; }
}

# load a sequence_set row, unless there is a set_id
if (! $set_id) {
    my $dbh = &connect({'u' => $user, 'p' => $pswd, 'D' => $db});
    $set_id = load_sequence_sets($dbh, $file_prefix, $organism, 0);
}

#cwd("/projects/ENVdbs/$db/");
#mkdir($shortName);
#cwd($shortName);

# load the genome sequence
if (! -e "${genome_fasta}.LOADED") {
    my $load_fasta_err = system("$ENV{ENVSCRIPTS}/load_fasta_to_ENVDB.pl -f $genome_fasta -D $db -p $pswd -s $set_id");
    if ($load_fasta_err) { die "Bad load_fasta_to_ENVDB.pl command\n" }
    else { system("touch ${genome_fasta}.LOADED") }
}


################################################################################
# Run prokka
################################################################################
my $prokka_gff = $file_prefix . ".gff";
if (! -e "${prokka_gff}.LOADED") {
    my $prokka_fsa = $file_prefix . ".fasta";
    my $prokka_faa = $file_prefix . ".faa";

    if (!( -e $prokka_gff &&
	   -e $prokka_fsa &&
	   -e $prokka_faa)) {

	my $prokka_gram = "";
	if ($gram) {
	    $prokka_gram = "--gram $gram";
	}
	my $prokka_org = "";
	if ($organism) {
	    my($genus, $species, $strain) = split(/\s+/,$organism,3);
	    $prokka_org = "--genus $genus";
	    $prokka_org .= " --species $species" if ($species);
	    $prokka_org .= " --strain $strain" if ($strain);
	}

	my $prokka_cmd = "$PROKKA --prefix $file_prefix --locustag $locus_prefix $prokka_org $prokka_gram --cpus $cpus --force $genome_fasta";
	print "$prokka_cmd\n";
	my $prokka_err = system($prokka_cmd);
	if ($prokka_err) { die "Bad Prokka command\n";}
    }

    if (! $protein_faa && -e $prokka_faa) { $protein_faa = $prokka_faa }
    elsif (! $protein_faa || -z $protein_faa) { die "Need to provide a protein fasta\n" }
    # Load prokka
    my $load_prokka_err = system("$ENV{ENVSCRIPTS}/gff_to_envdb.pl -D $db -u $user -p $pswd -g $prokka_gff -f $prokka_fsa -s PROKKA");
    if ($load_prokka_err) { die "Bad gff_to_envdb.pl command\n" }
    else { system("touch ${prokka_gff}.LOADED") }

}

if (! $protein_faa) {
    $protein_faa = $file_prefix . ".faa";
    system("$ENV{ENVSCRIPTS}/set_to_pep.pl -D $db -i $set_id -o $protein_faa");
    if (! -r $protein_faa || -S $protein_faa) {
	die("Problem with building protein set from $db/set_id: $set_id\n");
    }
}

if (! -e "${file_prefix}.TIGRfam.domtblout.LOADED") {
    # Run TIGRfams
    my $TIGRfam_err = system("$HMMSEARCH -o $file_prefix.TIGRfam.out --tblout $file_prefix.TIGRfam.tblout --domtblout $file_prefix.TIGRfam.domtblout --noali --acc --cut_nc --cpu $cpus $TIGRFAMDB $protein_faa");
    if ($TIGRfam_err) { die "Bad TIGRfam command\n" }

    # Load TIGRfams
    my $load_TIGRfam_err = system("$ENV{ENVSCRIPTS}/load_hmms_to_ENV.pl -D $db -U $user -P $pswd -i $file_prefix.TIGRfam.domtblout");
    if ($load_TIGRfam_err) { die "Bad load TIGRfams command\n" }
    else { system("touch ${file_prefix}.TIGRfam.domtblout.LOADED") }

    my $ann_by_hmm_err = system("$ENV{ENVSCRIPTS}/ann_by_hmm.pl -D $db -u $user -p $pswd -s $set_id");
}

if (! -e "${file_prefix}.Pfam.domtblout.LOADED") {
    # Run Pfam
    my $Pfam_err = system("$HMMSEARCH -o $file_prefix.Pfam.out --tblout $file_prefix.Pfam.tblout --domtblout $file_prefix.Pfam.domtblout --noali --acc --cut_ga --cpu $cpus $PFAMDB $protein_faa");
    if ($Pfam_err) { die "Bad Pfam command\n" }

    # Load Pfam
    my $load_Pfam_err = system("$ENV{ENVSCRIPTS}/load_hmms_to_ENV.pl -D $db -U $user -P $pswd -i $file_prefix.Pfam.domtblout");
    if ($load_Pfam_err) { die "Bad load pfams command\n" }
    else { system("touch ${file_prefix}.Pfam.domtblout.LOADED") }
}

if (! -e "${file_prefix}.KOfamscan.out.LOADED") {
    # run KOfams
    my $KOfamscan_err = system("$HMMSEARCH -o $file_prefix.KOfam.out --tblout $file_prefix.KOfam.tblout --domtblout $file_prefix.KOfam.domtblout --noali --acc -E 1e-25 --cpu $cpus $KOFAMDB $protein_faa");
    if ($KOfamscan_err) { die "Bad KOfamscan command\n" }

    # load KOfamscan
    my $load_KO_err = system("$ENV{ENVSCRIPTS}/load_KOfamscan.pl -D $db -u $user -p $pswd -i $file_prefix.KOfam.domtblout");
    if ($load_KO_err) { die "Bad load_KOfamscan command\n" }
    else { system("touch ${file_prefix}.KOfamscan.out.LOADED") }
}

if (! -e "./dbCAN.LOADED") {
    # Run dbCAN
    my $dbCAN_gram = $gram eq "pos" ? "p" : "n";
    my $dbCAN_err = system("$DBCAN $protein_faa protein --db_dir $DBCANDB --dia_cpu $cpus --hmm_cpu $cpus --hotpep_cpu $cpus --tf_cpu $cpus --stp_cpu $cpus --use_signalP T --gram $dbCAN_gram --out_dir dbCAN_out");
    if ($dbCAN_err) { die "Bad dbCAN command\n" }

    # Load dbCAN
    my $load_dbCan_err = system("$ENV{ENVSCRIPTS}/load_dbCAN_to_ENV.pl -D $db -u $user -p $pswd -d ./dbCAN_out");
    if ($load_dbCan_err) { die "Bad load_dbCAN_to_ENV.pl command\n" }
    else { system("touch dbCAN.LOADED") }
}

if (! -e "signalp.LOADED" && -e "./dbCAN_out/signalp.out") {
    my $load_signalp_err = system("$ENV{ENVSCRIPTS}/load_signalp.pl -D $db -u $user -p $pswd -i ./dbCAN_out/signalp.out");
    if ($load_signalp_err) { die "Bad load_signalp.pl command\n" }
    else { system("touch ${file_prefix}.signalp.LOADED") }
}

if (! -e "${file_prefix}.tmhmm.LOADED") {
    # Run tmhmm
    my $tmhmm_err = system("$TMHMM $protein_faa > $file_prefix.tmhmm");
    if ($tmhmm_err) { die "Bad tmhmm command\n" }

    # Load tmhmm
    my $load_tmhmm_err = system("$ENV{ENVSCRIPTS}/load_tmhmm.pl -D $db -u $user -p $pswd -i $file_prefix.tmhmm");
    if ($load_tmhmm_err) { die "Bad load_tmhmm command\n" }
    else { system("touch ${file_prefix}.tmhmm.LOADED") }
}

if (! -e "${file_prefix}.psortb.LOADED") {
    # Run psortb
    my $psort_gram = "";
    if ($gram) {
	$psort_gram = $gram eq "neg" ? "-n" : "-p";
    }

    my $psort_arch = "";
    if ($arch) {
	$psort_arch = "-a";
    }

    my $psort_err = system("$PSORT $psort_gram $psort_arch -o terse $protein_faa > $file_prefix.psortb");
    if ($psort_err) { die "Bad psort command\n" }

    # Load psortb
    my $load_psort_err = system("$ENV{ENVSCRIPTS}/load_psort.pl -D $db -u $user -p $pswd -i $file_prefix.psortb");
    if ($load_psort_err) { die "Bad load_psort.pl command\n" }
    else { system("touch ${file_prefix}.psortb.LOADED") }
}
