Winjao
======

Perl script for bioinformatics




NAME
       genbank_to_tbl - Convert a GBK file to Sequin Table format.

SYNOPSIS

       genbank2tbl.pl <-i/--gbk=INPUT>    <-o/--organism> [-s/--strain] [-l/--locus_tag] [-p/--prefix] [-h/--help].

DESCRIPTION

       This scripts generates files for ncbi sequence submition.  Output a Sequin Table format file (submit.tbl) and a Fasta format file (submit.fsa).

       example
       
           chaos_plot -i MYFIL.gbk -o "Shigella flexneri" -s "2a" -p "contig"


         You must give a specific organism name of the sequence.  If no strain name is given,


AUTHOR - Wenchao Lin
       Email linwenchao@yeah.net
