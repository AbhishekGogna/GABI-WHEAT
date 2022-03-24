#!/bin/bash
#$1=project directory # tells where the whole project is
#$2=name blastdb # tells where the database you made in that project is
#$3=query # tells where in the project the fasta file to be blasted is
#$4=output file # tells where to store the output in the project directory

ncbi-blast-2.12.0+/bin/blastn -num_threads 20 -db $1/$2 -query $1/$3 -evalue 1e-1 -dust no -out $1/$4_raw -outfmt '7 qseqid qlen qstart qend sseqid slen sstart send length pident mismatch gaps evalue bitscore sstrand'
echo "blastn done!"
grep -v -P "#" $1/$4_raw > $1/$4
echo "file formatted"