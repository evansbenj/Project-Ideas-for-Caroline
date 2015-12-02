# Developing a pipeline for analysis of abyss output

Or go back to the [previous page] (https://github.com/evansbenj/Project-Ideas-for-Caroline/blob/master/Project_9_rat_balls.md).

Caroline and Ben discussed a strategy for making sense of the abyss de novo assembly.  The plan is to take the de novo assembly (filename: ABTC26654-8.fa), which is a fasta file, and blast it against the mouse and rat genomes.  We will save the output as the top hit for each query sequence for each database.  

We think we will start by being 'liberal' with our designation of whether or not a region is on the X chromosome because we can easily test this later by looking for heterozygous sites in males, possibly from the shotgun sequences if the individual is male, which it is not for ABTC26654 but *probably* is for the other two shotgun sequences (check this, and definitely when we map the RADseq data which includes at least one male per species (I think).

We downloaded to the new computer NCBI Blast version 2.2.31+ and we downloaded the mouse and rat genomes we were using (mm10 and rn6 respectively) from sharcnet.

First we formatted the genome files like this:

`/Users/evanslab/ncbi-blast-2.2.31+/bin/makeblastdb -in /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked.fasta -dbtype nucl -title UCSC_mouse -out /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked_blastable`

for mouse, and like this:

`/Users/evanslab/ncbi-blast-2.2.31+/bin/makeblastdb -in /Users/Ben/rat_genomes/UCSC_rat_genome_rn6/rn6.masked.fa -dbtype nucl -title rat_rn6 -out /Users/Ben/rat_genomes/UCSC_rat_genome_rn6/rn6.masked_blastable`

for rat

And then we set up a batch blast to each one line this:

` /Users/evanslab/ncbi-blast-2.2.31+/bin/blastn -evalue 1e-60 -query /Users/Ben/rat_genomes/ABTC26654_abyss_kmer_65/ABTC26654-8.fa -db /Users/Ben/rat_genomes/UCSC_mouse_genome/mouse_genome_masked_blastable -out /Users/Ben/rat_genomes/ABTC26654_abyss_kmer_65/ABTC26654-8_mouse -outfmt 6 -max_target_seqs 1`

# Demultiplexing and trimming the RADseq data

Here is the data for the RADseq samples:
```
species	sample (RADseq)	Sex
Notomys	ABTC26646	Female
Notomys	ABTC26647	Female
Notomys	ABTC26648	Female
Notomys	ABTC26649	Female
Notomys	ABTC26652	Male
Notomys	ABTC26653	Male
Paruromys dominator 	JAE4305	Male
Paruromys dominator 	JAE4306	Male
Paruromys dominator 	JAE4307	Male
Paruromys dominator 	JAE4358	Female
Paruromys dominator 	JAE4399	Female
Paruromys dominator 	KCR1615	Female
Apodemus semotus	MVZ180233	Female
Apodemus semotus	MVZ180235	Female
Apodemus semotus	MVZ180237	Female
Apodemus semotus	MVZ180239	Female
Apodemus semotus	MVZ180241	Male
Apodemus semotus	MVZ180314	Male
Apodemus sylvaticus	SMG3824	Female
Apodemus sylvaticus	SMG3893	Male
Apodemus sylvaticus	SMG3894	Female
Apodemus sylvaticus	SMG902	Female
Apodemus sylvaticus	SMG3909	Female
Apodemus sylvaticus	SMG3910	Female
```

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

I then used Trimmomatic version 0.33 on to cleanup all samples using this perl script (Step_1_trimmomatic.pl):

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

Then I indexed the abyss assembly line this (Step_2_index_genome.pl):

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will index a genome fasta file using the 
# old bwa commands (not the HTSlib - this ended up stalling for GATK)

my $path_to_reference_genome="/home/ben/2015_rat_RADtags/reference_genomez_from_abyss/";
my $reference_genome="ABTC26654-8.fa";
my $path_to_picard = "/usr/local/picard-tools-1.131/";
my $status;
my $commandline;

# index the reference genome
$commandline = "bwa index  -a bwtsw ".$path_to_reference_genome.$reference_genome;
$status = system($commandline);

# make a fai index file for GATK
$commandline = "samtools faidx ".$path_to_reference_genome.$reference_genome;
$status = system($commandline);
# make a dict file for this genome
$commandline = "java -jar ".$path_to_picard."picard.jar CreateSequenceDictionary REFERENCE=".$path_to_reference_genome.$reference_genome." OUTPUT=".$path_to_reference_genome.$reference_genome.".dict";
$status = system($commandline);

# fix name of dict file
my $y = substr($reference_genome, 0, rindex($reference_genome, '.'));
$commandline="mv ".$path_to_reference_genome.$reference_genome.".dict ".$path_to_reference_genome.$y.".dict";
$status = system($commandline);
````
And then I mapped the reads like this (Step_3_align_RADtags.pl):

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will index a genome fasta file using the 
# new bwa HTSlib commands  

my $path_to_reference_genome="/home/ben/2015_rat_RADtags/reference_genomez_from_abyss/";
my $reference_genome="ABTC26654-8.fa";
my $path_to_GATK="/usr/local/gatk/";
my $path_to_picard = "/usr/local/picard-tools-1.131/";
my $commandline;
my $status;

my @files = glob("ABTC*fastq.gz");

#foreach(@files){
#    # align the data for each individual
#    $commandline = "bwa mem -M -t 4 ".$path_to_reference_genome.$reference_genome." -R \'\@RG\\tID:".$_."\\tSM:".$_."\\tLB:library1\\tPL:illumina\' ".$_." | gzip -3 > ".$_.".sam";
#    $status = system($commandline);
#    # convert the sam file to a sorted bam file with picard
#    $commandline = "java -jar ".$path_to_picard."picard.jar SortSam INPUT=".$_.".sam OUTPUT=".$_.".sorted.bam SORT_ORDER=coordinate";
#    $status = system($commandline);
#    # delete the sam file
#    $commandline = "rm -f ".$_.".sam";
#    $status = system($commandline);
#    $commandline = "java -jar ".$path_to_picard."picard.jar BuildBamIndex INPUT=".$_.".sorted.bam";                       
#    $status = system($commandline);                                                            
                           
#    $commandline = "cp ".$_.".sorted.bam.bai ".$_."sorted.bai";
#    $status = system($commandline);
#}

# Now identify Indels
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T RealignerTargetCreator
 -nt 8";
foreach(@files){
    $commandline = $commandline." -I ".$_."sorted.bam ";
}
$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." -o indel.intervals";
print $commandline,"\n";
$status = system($commandline);

# Now realign indels
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T IndelRealigner ";
foreach(@files){
    $commandline = $commandline." -I ".$_."sorted.bam ";
}
$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." --targetIntervals indel.intervals --nWayOut.realigned.bam";
$status = system($commandline);

# Now recalibrate bases; first emit nonrecal variants
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline." -out_mode EMIT_VARIANTS_ONLY -o ./nonrecal_varonly.vcf";
$status = system($commandline);

# Now do baserecalibration
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T BaseRecalibrator -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline." -knownSites ./nonrecal_varonly.vcf -o recal.table";
$status = system($commandline);

# Print new concatenated recalibrated bam
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T PrintReads -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline."-BQSR recal.table -o all_recal_round1.bam";
$status = system($commandline);

# Now call all sites
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." -I all_recal_round1.bam ";
$commandline = $commandline." -out_mode EMIT_ALL_CONFIDENT_SITES -o ./recal_allsites.vcf";
$status = system($commandline);

```

The problem is that with the trial using the ABCT genome assembly, the pipeline is stalling at the "RealignerTargetCreator" stage. Not sure why.
