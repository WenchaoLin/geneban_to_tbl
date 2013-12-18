#!/usr/bin/perl

use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqFeatureI;
use Data::Dumper;

my $gbk;
my $organism;
my $strain;
my $output = "submit";
my $locus_tag="LOC";
my $topology="liner";
my $pre="ctg";
my $help=0;
my $gcode = 11;

my $description = "[molecule=DNA] [tech=wgs]";


GetOptions ("gbk|i:s" => \$gbk,
                "organism|o:s" =>\$organism,
                "strain|S:s" =>\$strain,
                "output:s" =>\$output,
                "gcode:i" =>\$gcode,
                "topology:s" =>\$topology,
                "locus_tag|l:s" =>\$locus_tag,
                "prefix|p:s" =>\$pre,
                "help|h" =>\$help,
           );

$USAGE = "usage:\tgenbank2tbl.pl <-i/--gbk=INPUT> <-o/--organism>\n\n".
"\t<-i/--gbk> \tinput genbank format file\n".
"\t<-o/--organism> \torganism name\n".
"\t-p/--output prefix \toutput prefix for fsa and tbl. Default:submit\n".
"\t-s/--strain \tstrain name\n".
"\t--topology \t[topology=liner|circular] (default is liner)\n".
"\t--gcode \tgenetic code. [gcode=11] (default is 11)\n".
"\t-l/--locus_tag \tlocus_tag you want to use. (default is LOC)\n".
"\t-p/--prefix \t[default:ctg] Contig name prefix\n".
"\t-h/--help\tshow this message.\n";

$USAGE .= "\tConvert a GBK file to Sequin Table format.\n";
$USAGE .= "\tDefault ouput will be in submit.tbl and submit.fsa\n";
$USAGE .= "\tBug report to Wenchao Lin <linwenchao\@yeah.net>\n\n";


if( ! $gbk || $help || !$organism) {
        die $USAGE ;
}


if (defined $organism){
        $description .= " [organism=$organism]";
}
if (defined $strain){
        $description .= " [strain=$strain]";
}
if (defined $gcode){
        $description .= " [gcode=$gcode]";
}
if (defined $topology){
        $description .= " [topology=$topology]";
}


my @allowed_tags = ('locus_tag',  'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag');
#my @allowed_tags = ('locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag');


open TBL, ">$output.tbl" or die "open $output.tbl error $!\n";
my $seqin = Bio::SeqIO->new(-file =>$gbk,-format=>'genbank');
my $seqout = Bio::SeqIO->new(-file=>">$output.fsa",-format=>'fasta');

my $count = 0; 
my $last;
my $tag_count = 0;

while(my $seq=$seqin->next_seq()){
        my $seq_id = $seq->display_id();
        my $seq_len = $seq->length();

        if ($seq_len <=200)
        {
                print STDERR "$seq_id length is $seq_len. less than 200, skipping small contig...\n";
                next;
        }

        $count += 1;
        my $id = $pre . sprintf ("%04d",$count);

        $seq->seq($seq->seq); 
        $seq->display_id($id);
        $seq->desc($description);

        $seqout->write_seq($seq);

        printf TBL ">Feature %s\n",$id;

        my %hash;

        my @feats = $seq->get_all_SeqFeatures();
        for my $f (@feats){

                my $tag = $f->primary_tag();
                my $start = $f->start();
                my $end = $f->end();
                my $strand = $f->strand();

                my @tmp = $f->get_tag_values('locus_tag');
                if ($tmp[0] ne $last){
                        $last = $tmp[0];
                        $tag_count += 1;
                }


                my $new_tag = $locus_tag . sprintf ("%04d",$tag_count);

                foreach my $key (@allowed_tags){
                        delete $hash{$key};
                        if($f->has_tag($key)){
                                my @value = $f->get_tag_values($key);
                                if($key eq "locus_tag"){
                                        $hash{$key} = $new_tag; 
                                }else{
                                        $hash{$key} = $value[0];
                                }
                        }
                }


                if($strand == 1){
                        printf TBL ("%d\t%d\t%s\n",$start,$end,$tag);
                }
                else
                {
                        printf TBL ("%d\t%d\t%s\n",$end,$start,$tag);
                }

                foreach (keys %hash){
                        printf TBL "\t\t\t%s\t%s\n",$_,$hash{$_};
                }

        }
}


__END__

=head1 NAME

genbank_to_tbl - Convert a GBK file to Sequin Table format.

=head1 SYNOPSIS

genbank2tbl.pl <-i/--gbk=INPUT>	<-o/--organism> [-s/--strain] [-l/--locus_tag] [-p/--prefix] [-h/--help].

=head1 DESCRIPTION

This scripts generates files for ncbi sequence submition.
Output a Sequin Table format file (submit.tbl) and a Fasta
format file (submit.fsa).

=over 2

=item example

  chaos_plot -i MYFIL.gbk -o "Shigella flexneri" -s "2a" -p "contig"

==item notice

You must give a specific organism name of the sequence.
If no strain name is given, 

=back

=head2 Reporting Bugs

Bug reports can be submitted via the email:

  linwenchao@yeah.net

=head1 AUTHOR - Wenchao Lin

Email linwenchao@yeah.net

=head1 HISTORY

Written by Wenchao Lin (TBC) linwenchao@yeah.net

=cut


