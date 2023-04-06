#!/bin/bash

# input is the genome multi-fasta file
# output will be a file of protein sequences from each genome
infasta=$1
filebase=`basename $infasta`
cpu=8

# infeRNAl/cmscan/Rfam:
# infernal can run on many file formats automagically, and will use as many processors as are available.
cmscanout=${filebase/.fna/.cmscan}
cmscantblout=$cmscanout.tblout
clanin=/projects/db/hmms/Rfam/Rfam.clanin
Rfamdb=/projects/db/hmms/Rfam/Rfam.cm

if [ ! -s "$cmscantblout" ]; then
    echo "running cmscan on $infasta"
    cmscan --noali --notextw --cut_ga --rfam --nohmmonly --cpu $cpu --tblout $cmscantblout --fmt 2 --clanin $clanin $Rfamdb $infasta > $cmscanout
fi

# Output will include models we're not interested in, so make a list of types that should be masked.
# Includes rRNAs and tRNAs. Rfam doesn't identify which tRNA, but sufficient for masking.
# From infernal output, parse the rRNA and tRNA (unless the anyidx field is populated)
# and mask the genome file.
#Generates the .masked file.
maskedfasta=${filebase/.fna/_masked.fna}

if [ ! -s "$maskedfasta" ]; then
    echo "masking $infasta"
    cat $cmscantblout | perl -ne 'if (!/^#/) { @f=split(/\s+/,$_, 27); if ($f[1] =~ /(t|r)RNA/) { print "$f[3]\t$f[9]\t$f[10]\n" unless ($f[20] =~ /\d+/);}}' | /projects/perlscripts/mask_fasta.pl $infasta
fi

# Okay now for the gene callers.
# Prodigal is easy/fast
PRODIGAL=${maskedfasta/.fna/.Prodigal.gff}

if [ ! -s "$PRODIGAL" ]; then
    echo "running Prodigal on $maskedfasta"
    prodigal -f gff -g 11 -i $maskedfasta -m -o $PRODIGAL -p single 
fi

# Glimmer3
GLIMMER=${maskedfasta/.fna/.Glimmer3}
if [ ! -s "$GLIMMER.predict" ]; then
    echo "running glimmer on $maskedfasta"
    /projects/bifx/glimmer3.02/scripts/multifasta_g3_iterated.py --infile $maskedfasta --tag $GLIMMER
fi

#GenemarkS - actually seems to use the same basic calculation as Prodigal, but seems to be more stringent in gene calls
GENEMARK=${maskedfasta/.fna/.GeneMarkS.gff}
if [ ! -s "$GENEMARK" ]; then
    echo "running geneMark on $maskedfasta"
    /projects/bifx/gms2_linux_64/gms2.pl --seq $maskedfasta --genome-type bacteria --format gff3 --output $GENEMARK
fi

#Now integrate the three gene callers with integrate_gene_calls.py
genecalls=${filebase/.fna/.CDS.tsv}
if [ ! -s "../proteins/${filebase/.fna/.protein.faa}" ]; then
    python3 /projects/pythontools/integrate_gene_calls.py --prodigal $PRODIGAL --glimmer $GLIMMER.predict --genemark $GENEMARK --fasta $maskedfasta --output ../gene_finding/$genecalls --gene-fasta ../genes/${filebase/.fna/.gene.fna} --protein-fasta ../proteins/${filebase/.fna/.protein.faa}
fi


