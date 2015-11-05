# Developing a pipeline for analysis of abyss output

Or go back to the [previous page] (https://github.com/evansbenj/Project-Ideas-for-Caroline/blob/master/Project_9_rat_balls.md).

Caroline and Ben discussed a strategy for making sense of the abyss de novo assembly.  The plan is to take the de novo assembly (filename: ABTC26654-8.fa), which is a fasta file, and blast it against the mouse and rat genomes.  We will save the output as the top hit for each query sequence for each database.  

We think we will start by being 'liberal' with our designation of whether or not a region is on the X chromosome because we can easily test this later by looking for heterozygous sites in males, possibly from the shotgun sequences if the individual is male, which it is not for ABTC26654 but *probably* is for the other two shotgun sequences (check this, and definitely when we map the RADseq data which includes at least one male per species (I think).

We downloaded to the new computer NCBI Blast version 2.2.31+ and we downloaded the mouse and rat genomes we were using (XXX and rn6 respectively) from sharcnet.

First we formatted the genome files like this:

`/Users/evanslab/ncbi-blast-2.2.31+/bin/makeblastdb -in /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked.fasta -dbtype nucl -title UCSC_mouse -out /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked_blastable`

for mouse, and like this:

`/Users/evanslab/ncbi-blast-2.2.31+/bin/makeblastdb -in /Users/Ben/rat_genomes/UCSC_rat_genome_rn6/rn6.masked.fa -dbtype nucl -title rat_rn6 -out /Users/Ben/rat_genomes/UCSC_rat_genome_rn6/rn6.masked_blastable`

for rat

And then we set up a batch blast to each one line this:

` /Users/evanslab/ncbi-blast-2.2.31+/bin/blastn -evalue 1e-60 -query /Users/Ben/rat_genomes/ABTC26654_abyss_kmer_65/ABTC26654-8.fa -db /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked_blastable -out /Users/Ben/rat_genomes/ABTC26654_abyss_kmer_65/ABTC26654-8_mouse -outfmt 6 -max_target_seqs 1`

# Demultiplexing and trimming the RADseq data

We initially demultiplexed the RADseq data using RADpools version 1.2.4.  
Here is the command we used:
`RADpools -i C213_1.fastq -d rat -s -f -v -o`

This uses the fuzzy barcode option (-f) which assigns missing barcodes to the nearest pool (or an unassigned pool for ties), outputs fastq files (-o), defines the quality scores as sanger (-s), and does a verbose output (-v).

Here is the barcode file (referenced by -d rat flag above):

```
MVZ180233 TCCGGAGCGC
MVZ180235 CTAACACGGC
MVZ180237 AGCTTCGATT
MVZ180239 TCGCCGCAAT
MVZ180241 TCAGTTCCGG
MVZ180314 CGGAAGTGAG
ABTC26646 GTTGCTAGAC
ABTC26647 AATAGATTCA
ABTC26648 AGCTGATACA
ABTC26649 ATCAGTAGAA
ABTC26652 TCGTCTTAGT
ABTC26653 GCTCAGCCAG
SMG3824 CGGCTACTTC
SMG3893 CAAGCCGGTT
SMG3894 TTGCGCAAGC
SMG3902 TACGATGGAG
SMG3909 GCAATATACA
SMG3910 AAGAATTCGG
KCR1615 TCGGCAGTCG
JAE4307 AGTTCCATTG
JAE4358 TTCTTGCGCT
JAE4399 AGCAATCTAA
JAE4305 GAATTGTCGC
JAE4306 CTTCGACATA
```

Using fastQC version 0.11.3, I identified these repetitive sequences in the data (pooled across all samples):
```
>repeat1
AGCAATCTAATGCAGGAGCTCCAGCATGGATAGGAGAGGGGCTCCCACAC
>repeat2
AGCAATCTAATGCAGGTCTTGGACAGGAAACCACAGCTGCTGTGTAATAT
>repeat3
CTTCGACATATGCAGGAGCTCCAGCATGGATAGGAGAGGGGCTCCCACAC
>repeat4
CTTCGACATATGCAGGTCTTGGACAGGAAACCACAGCTGCTGTGTAATAT
>repeat5
GAATTGTCGCTGCAGGAGCTCCAGCATGGATAGGAGAGGGGCTCCCACAC
>repeat6
GAATTGTCGCTGCAGGTCTTGGACAGGAAACCACAGCTGCTGTGTAATAT
>repeat7
TGCAGGAAAGGCAAGACACAACCAGACAGAGTAAGATCAACTAACATCAG
>repeat8
TGCAGGAACACATCTTTGCAGACCTAGTTTGTGGCTACGCACAACTTCTC
>repeat9
TGCAGGAACACATCTTTGCAGACCTAGTTTGTGGCTACGCACAACTTTTC
>repeat10
TGCAGGAAGTAACTGAGGCCCCAGCAGCTCCCCTTCCACTTAAAGCAAGT
>repeat11
TGCAGGAATACTAGAAGTCCTGACTATGCCACCCAACTGGATGGCCAACA
>repeat12
TGCAGGACCTCCTGGTACTAGGAGACAAGGGTTCCCCATCATGACCCAAG
>repeat13
TGCAGGACGACAGAGAGCTCTGGCGGGCAGAGCCATAACAGGCATAGCCT
>repeat14
TGCAGGACTACAGAGAGCTCGGGCGGGCAGAGCCATAACAGGCATAGCCT
>repeat15
TGCAGGAGACTCAGGCAAGCTGCACACCCCCACATCCTGCAGCCACCACT
>repeat16
TGCAGGAGACTCTGGCAAGCTGCACACCCCCACATCCTGCAGCCACCACT
>repeat17
TGCAGGAGCTCCAGCATGGAGAGGAGAGGGGCTCCCACACCTGTGGAAGA
>repeat18
TGCAGGAGCTCCAGCATGGATAGGAGAGGGGCTCCCACACCTGTGGAAGA
>repeat19
TGCAGGAGCTCCAGCATGGATAGGAGAGGGGCTCCTACACCTGTGGAAGA
>repeat20
TGCAGGAGGGACAAGACACAACCAGACAGAGTAAGATCAACTAACATCTG
>repeat21
TGCAGGAGGGACAAGACACAACCAGACAGAGTAAGATCAACTAACTTCAG
>repeat22
TGCAGGAGGGACAAGACACATCCAGACAGAGTAAGATCAACTAACATCTG
>repeat23
TGCAGGAGGGACAAGCCACAGCCAGACAGAGTAAGATCAGCTAACATCAG
>repeat24
TGCAGGAGGGGCAAGACACAACCAGACAGAGTAAGATCAACTAACATCAG
>repeat25
TGCAGGAGGGGCAAGCCACAACCAGACAGAGTAAGATCAGCTAACATCAG
>repeat26
TGCAGGCCAGAAGATTATTATAGACGTTGTTCCTTATCTCATTACAGATG
>repeat27
TGCAGGCCTGGAACTTGCTCTGTAGACTTCATTAGCCTCAAACTCACAGT
```

I then used Trimmomatic version 0.33 on to cleanup all samples using this perl script:

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will run trimmomatic on SE radtag reads   

my $status;
my $commandline;
@files = glob("*.fastq");

foreach(@files){
	$commandline = "java -jar /usr/local/trimmomatic/trimmomatic-0.33.jar SE -trimlog ".$_."txt ".$_." ".$_."single.fastq.gz ILLUMINACLIP:/home/ben/2015_rat_RADtags/adapters/RADseq_adaptors.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36";
	$status = system($commandline);
}
```
