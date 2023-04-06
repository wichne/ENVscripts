# ENVscripts
Perl scripts that interact with the ENV relational database structure, loading, retrieving, etc

To load a genome:
1) First make a sequence_set entry
$ for f in *.gtf; do g=${f/\.gtf/}; echo "insert into sequence_sets (name, description, is_current) values (\"$g\", \"Neisseria $g\", 1)"; done | runsql -D Neisseria -P RenMan

2) Load the genome fasta file.
