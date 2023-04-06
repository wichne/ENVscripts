#!bash
while read set; do
    perl /share/scripts/devel/ENV/set_to_fasta.pl -D OD1 -n $set -o $set.fsa
    tRNAscan-SE $set.fsa > $set.fsa.tRNAscan
    /share/scripts/devel/ENV/load_tRNAs_to_ENV.pl -d OD1 -p RenMan -i $set.fsa.tRNAscan
    echo "Done with $set tRNAs."
    rfam_scan.pl -o $set.fsa.rfam -blastdb /common/db/hmms/Rfam/Rfam.fasta --nobig --exclude RF00005 --exclude RF01852 /common/db/hmms/Rfam/Rfam.cm $set.fsa
    perl /share/scripts/devel/ENV/load_Rfam_to_ENV.pl -d OD1 -p RenMan -i $set.fsa.rfam 
done;
