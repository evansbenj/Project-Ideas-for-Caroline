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
Notomys	mitchelli ABTC26646	Female
Notomys	mitchelli ABTC26647	Female
Notomys	mitchelli ABTC26648	Female
Notomys	mitchelli ABTC26649	Female
Notomys	mitchelli ABTC26652	Male
Notomys	mitchelli ABTC26653	Male
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
Apodemus sylvaticus	SMG3902	Female
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


Based on a brief comment in the GATK forum, it is now clear to me that GATK was not designed to work with reference genomes with many, many contigs.  So we now have a quick fix.  Based on the BLAST results of the abyss assembly to the mouse and rat genomes, we are making "supercontigs" using this perl script (Combines_abyss_output_into_supercontigs.pl). Note that GATK also cannot handle really really large supercontigs, so I have split up the contig that map to autosomes in mouse and rat into multiple contigs, each with ~1 million lines of data:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


# This program will read in a several lists of contig names
# that were generated based on a blast run against mouse and
# rat genomes.  It will then make a master fastafile that includes
# concatenated reads from the following types of contigs:
# (1) contigs with a top blast hit to mouse chrX and rat chrX
# (2) contigs with a top blast hit to mouse chrX but a rat autosome
# (3) contigs with a top blast hit to rat chrX but a mouse autosome
# (4) contigs with a top blast hit to mouse chrX and a rat unknown region
# (5) contigs with a top blast hit to mouse autosome and a rat autosome
# (6) contigs with a top blast hit to mouse autosome and a rat unknown region

# Combines_abyss_output_into_supercontigs.pl ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Xrat_ABTC.txt ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Arat_ABTC.txt ./abyss_contig_blast_lists/ABTC/Amouse_AND_Xrat_ABTC.txt ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Urat_ABTC.txt ./abyss_contig_blast_lists/ABTC/Amouse_AND_Arat_ABTC.txt ./abyss_contig_blast_lists/ABTC/Amouse_AND_Urat_ABTC.txt ../reference_genomez_from_abyss/ABTC26654-8.fa ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/ABTC/Amouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/ABTC/Xmouse_AND_Urat_ABTC.bed ./abyss_contig_blast_lists/ABTC/Amouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/ABTC/Amouse_AND_Urat_ABTC.bed ../reference_genomez_from_abyss/ABTC26654-8_concat.fa
# Combines_abyss_output_into_supercontigs.pl ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Xrat_MVZ.txt ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Arat_MVZ.txt ./abyss_contig_blast_lists/MVZ/Amouse_AND_Xrat_MVZ.txt ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Urat_MVZ.txt ./abyss_contig_blast_lists/MVZ/Amouse_AND_Arat_MVZ.txt ./abyss_contig_blast_lists/MVZ/Amouse_AND_Urat_MVZ.txt ../reference_genomez_from_abyss/MVZ180318-8.fa ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/MVZ/Amouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/MVZ/Xmouse_AND_Urat_ABTC.bed ./abyss_contig_blast_lists/MVZ/Amouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/MVZ/Amouse_AND_Urat_ABTC.bed ../reference_genomez_from_abyss/MVZ180318-8_concat.fa
# Combines_abyss_output_into_supercontigs.pl ./abyss_contig_blast_lists/JAE/Xmouse_AND_Xrat_JAE.txt ./abyss_contig_blast_lists/JAE/Xmouse_AND_Arat_JAE.txt ./abyss_contig_blast_lists/JAE/Amouse_AND_Xrat_JAE.txt ./abyss_contig_blast_lists/JAE/Xmouse_AND_Urat_JAE.txt ./abyss_contig_blast_lists/JAE/Amouse_AND_Arat_JAE.txt ./abyss_contig_blast_lists/JAE/Amouse_AND_Urat_JAE.txt ../reference_genomez_from_abyss/JAE4405-8.fa ./abyss_contig_blast_lists/JAE/Xmouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/JAE/Xmouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/JAE/Amouse_AND_Xrat_ABTC.bed ./abyss_contig_blast_lists/JAE/Xmouse_AND_Urat_ABTC.bed ./abyss_contig_blast_lists/JAE/Amouse_AND_Arat_ABTC.bed ./abyss_contig_blast_lists/JAE/Amouse_AND_Urat_ABTC.bed ../reference_genomez_from_abyss/JAE4405-8_concat.fa

my $minimum_length_of_contig_to_include = 500;
my $max_number_of_lines_in_autosomal_supercontigs=2000000;
my $achrom_lines=0;
my $a_chr_number=0;

my @input;
$input[0] = $ARGV[0]; #This is list 1 defined above >chrX_m_chrX_r
$input[1] = $ARGV[1]; #This is list 2 defined above >chrX_m_chrA_r
$input[2] = $ARGV[2]; #This is list 3 defined above >chrA_m_chrX_r
$input[3] = $ARGV[3]; #This is list 4 defined above >chrX_m_chrU_r
$input[4] = $ARGV[4]; #This is list 5 defined above >chrA_m_chrA_r
$input[5] = $ARGV[5]; #This is list 6 defined above >chrA_m_chrU_r
$input[6] = $ARGV[6]; #This is the input abyss assembly

my @output;
$output[0] = $ARGV[7]; # This is the bedfile for >chrX_m_chrX_r
$output[1] = $ARGV[8]; # This is the bedfile for >chrX_m_chrA_r
$output[2] = $ARGV[9]; # This is the bedfile for >chrA_m_chrX_r
$output[3] = $ARGV[10]; # This is the bedfile for >chrX_m_chrU_r
$output[4] = $ARGV[11]; # This is the bedfile for >chrA_m_chrA_r
$output[5] = $ARGV[12]; # This is the bedfile for >chrA_m_chrU_r

$output[6] = $ARGV[13]; # This is the fastaoutputfile 

my $y;
my $x;
my $counter;
my @temp;
my %list;
my @number_of_entries;

for ($y = 0 ; $y < 6 ; $y++ ) {
	$counter=0;
	unless (open DATAINPUT, $input[$y]) {
		print "Can not find the input file $y.\n";
		exit;
	}
	while ( my $line = <DATAINPUT>) {
		@temp=split(/\s+/,$line);
		$list{$y}[$temp[0]]="keep_me";
		$counter+=1;
	}
	$number_of_entries[$y]=$counter-1;
}



# open the concat fasta output file
unless (open(OUTFILE6, ">$output[6]"))  {
	print "I can\'t write to $output[6]\n";
	exit;
}
print "Creating output file: $output[6]\n";

# now, for each list, make the concatenated chromosome
# and print it to the outfile
# and make a bed file

my $switch=0;
my $last_position=0;
my $seq_length=0;
for ($y = 0 ; $y < 6 ; $y++ ) {
	$switch=0;
	$last_position=0;
	unless (open(OUTFILE, ">$output[$y]"))  {
		print "I can\'t write to $output[$y]\n";
		exit;
	}
	print "Creating output file: $output[$y]\n";
	
	if($y == 0){
		print OUTFILE6 ">chrX_m_chrX_r\n";
		$achrom_lines=0;
	}
	elsif($y == 1){
		print OUTFILE6 ">chrX_m_chrA_r\n";
		$achrom_lines=0;
	}
	elsif($y == 2){
		print OUTFILE6 ">chrA_m_chrX_r\n";
		$achrom_lines=0;
	}
	elsif($y == 3){
		print OUTFILE6 ">chrX_m_chrU_r\n";
		$achrom_lines=0;
	}
	elsif($y == 4){
		print OUTFILE6 ">chrA_m_chrA_r_".$a_chr_number."\n";
		$achrom_lines=0;
	}
	elsif($y == 5){
		print OUTFILE6 ">chrA_m_chrU_r\n";
		$achrom_lines=0;
	}

	# open the abyss assembly file
	unless (open DATAINPUT6, $input[6]) {
		print "Can not find the abyss assembly input file $y.\n";
		exit;
	}
	
	# now cycle through and print out the sequences from
	# only those contigs in the list for this chr	
	while ( my $line = <DATAINPUT6>) {
		@temp=split(/[>\s]/,$line);
		if(($switch == 1)&&($line !~ /^>/)){
			until(length($temp[0]) < 80){
				print OUTFILE6 substr($temp[0], 0, 80),"\n";
				$temp[0] = substr($temp[0],80);
				$achrom_lines+=1;
			}
			print OUTFILE6 $temp[0];
			for($x = 0 ; $x < (80-length($temp[0])) ; $x++ ) {
				print OUTFILE6 "N";
			}
			print OUTFILE6 "\n";
			for($x = 0 ; $x < 80 ; $x++ ) {
				print OUTFILE6 "N";
			}	
			print OUTFILE6 "\n";
			print OUTFILE $last_position+$seq_length+(160-length($temp[0])),"\n";
			$last_position = $last_position+$seq_length+(160-length($temp[0]));
			$switch=0;
		}
		elsif($line =~ /^>/){
			if((defined($list{$y}[$temp[1]]))&&($temp[2]>$minimum_length_of_contig_to_include)){
				$switch = 1;
				$seq_length=$temp[2];
				print OUTFILE $temp[1],"\t",$last_position+1,"\t";
				if(($achrom_lines > $max_number_of_lines_in_autosomal_supercontigs)&&($y == 4)){
					$a_chr_number+=1;
					$achrom_lines=0;
					$last_position=0;
					print OUTFILE6 ">chrA_m_chrA_r_".$a_chr_number."\n"; 
				}
			}
			else{
				$switch = 0;
			}
		}
	}
	close DATAINPUT6;
	#print OUTFILE6 "\n";
}



```

Then I indexed the supercontigs from abyss (Step_2_index_genome.pl):

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will index a genome fasta file using the 
# new bwa HTSlib commands  

my $path_to_reference_genome="/home/ben/2015_rat_RADtags/reference_genomez_from_abyss/";
my $reference_genome="MVZ180318-8_concat.fa";
my $path_to_picard = "/usr/local/picard-tools-1.131/";
my $status;
my $commandline;

# index the reference genome
$commandline = "bwa index ".$path_to_reference_genome.$reference_genome;
$status = system($commandline);

# make a fai index file for GATK
$commandline = "samtools faidx ".$path_to_reference_genome.$reference_genome;
$status = system($commandline);

# make a dict file for this genome
$commandline = "java -jar ".$path_to_picard."picard.jar CreateSequenceDictionary REFERENCE=".$path_to_reference_genome.$reference_genome." OUTPUT=".$path_to_reference_genome.$reference_genome.".dict";
print $commandline,"\n";
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
my $reference_genome="JAE4405-8_concat.fa";
my $path_to_GATK="/usr/local/gatk/";
my $path_to_picard = "/usr/local/picard-tools-1.131/";
my $path_to_bwa="/usr/local/bin/";
my $path_to_samtools="/usr/local/bin/";
my $commandline;
my $status;

my @files = glob("JAE*fastq.gz KCR*fastq.gz");



foreach(@files){
    # align the data for each individual
    $commandline = "bwa mem -M -t 4 ".$path_to_reference_genome.$reference_genome." -R \'\@RG\\tID:".$_."\\tSM:".$_."\\tLB:library1\\tPL:illumina\' ".$_." | gzip -3 > ".$_.".sam";
    print $commandline,"\n";
    $status = system($commandline);

    # convert the sam file to a sorted bam file with picard
    $commandline = "java -jar ".$path_to_picard."picard.jar SortSam INPUT=".$_.".sam OUTPUT=".$_.".sorted.bam SORT_ORDER=coordinate";
    print $commandline,"\n";
    $status = system($commandline);

    # delete the sam file
    $commandline = "rm -f ".$_.".sam";
    print $commandline,"\n";
    $status = system($commandline);

     # build an index
    $commandline = "java -jar ".$path_to_picard."picard.jar BuildBamIndex INPUT=".$_.".sorted.bam";    
    print $commandline,"\n";
    $status = system($commandline);                                                            
}



# Now identify Indels
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T RealignerTargetCreator ";
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.bam ";
}
$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." -o ".$reference_genome."_indel.intervals";
print $commandline,"\n";
$status = system($commandline);

# Now realign indels
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T IndelRealigner ";
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.bam ";
}
$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." --targetIntervals ".$reference_genome."_indel.intervals --nWayOut .realigned.bam";
print $commandline,"\n";
$status = system($commandline);


# Now recalibrate bases; first emit nonrecal variants
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline." -out_mode EMIT_VARIANTS_ONLY -o ".$reference_genome."_nonrecal_varonly.vcf";
print $commandline,"\n";
$status = system($commandline);

# Now do baserecalibration
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T BaseRecalibrator -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline." -knownSites ".$reference_genome."_nonrecal_varonly.vcf -o ".$reference_genome."_recal.table";
print $commandline,"\n";
$status = system($commandline);

# Print new concatenated recalibrated bam
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T PrintReads -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_.".sorted.realigned.bam ";
}
$commandline = $commandline."-BQSR ".$reference_genome."_recal.table -o ".$reference_genome."_all_recal_round1.bam";
print $commandline,"\n";  
$status = system($commandline);

# Now call all sites
$commandline = "java -Xmx3G -jar ".$path_to_GATK."GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." -I ".$reference_genome."_all_recal_round1.bam ";
$commandline = $commandline." -out_mode EMIT_ALL_CONFIDENT_SITES -o ".$reference_genome."_recal_allsites.vcf";
print $commandline,"\n"; 
$status = system($commandline);
```

After this, we are ready to convert the vcf file into a tab file using vcftools.  Here are the commands for that:

```bash
~/tabix-0.2.6/bgzip JAE4405-8_concat.fa_recal_allsites.vcf
~/tabix-0.2.6/tabix -p vcf JAE4405-8_concat.fa_recal_allsites.vcf.gz
zcat JAE4405-8_concat.fa_recal_allsites.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > JAE4405-8_concat.fa_recal_allsites.vcf.gz.tab
```

Then I wrote a script to identify heterozygous sites in males.  These and the contigs from which they originate will be deleted from the analysis (FInds_heterozygous_sites_in_supercontigs.pl):

```perl

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by vcftools and identifies sites that are heterozygous in males
#  and that are also on the X chromosome.  These sites will then be
#  used to identify and exclude entire contigs that are not actually
#  on the X.

# to execute type Finds_hetero.pl inputfile.tab 11100 output.tab
# where 11100 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female 

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];
my $y;
my @temp;


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my @sexes = split("",$ARGV[1]);

my $number_of_individuals_genotyped=($#sexes + 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped;
for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	

print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $switch=0;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		if(($temp[0] eq "chrX_m_chrX_r")||($temp[0] eq "chrX_m_chrA_r")||($temp[0] eq "chrA_m_chrX_r")||($temp[0] eq "chrX_m_chrU_r")){
			for ($y = 0 ; $y <= $#sexes; $y++ ) {
				# load both alleles if the individual is a male
				if($sexes[$y] eq "0"){
					# load the first allele
					if($temp[(2*$y)+3] ne $temp[(2*$y)+4]){
						$switch = 1;
					}	
				}
			} 
			if($switch == 1){
				print OUTFILE $temp[0],"\t",$temp[1],"\n";
			}
			$switch=0;
		}
	}
}	
close DATAINPUT;
close OUTFILE;
```


I sent a list of heterozygous sites for each of the 4 species to Caroline.  Caroline then wrote an R script that generates bed files with the coordinates of contigs that had heterozygous sites in the RADseq data in one or more males on the X chromosome.  These contigs are either not on the X or are pseudoautosomal and we need to eliminate them from our analysis. For the MVZ genome, we have a list of male X-heterozygous sites for semotus and sylvaticus that we will combine. We can then use vcf tools to delete these sections from our supercontigs.

To obtain a file with the "chromosome situation" in the fist column, the start and the end position of the contig with heterozygous sites (here is a script example used for ABTC):

``` R
## import_multiple_files_to_R
# set working directory
setwd("/Users/Ben/rat_genomes/ABTC")
#path = "C:/Users/Ben/rat_genomes/ABTC/"

###Import the bed files and the txt files containing the list of sites
filelist_bed = list.files(pattern = ".*.bed")

#read the heterozygous file
heteroz <- read.table("mit_X_het_sites.out",sep="\t")
head(heteroz)
levels(heteroz[,1])
chrA_m_chrX_r <- subset(heteroz,heteroz[,1] == "chrA_m_chrX_r")
head(chrA_m_chrX_r)
tail(chrA_m_chrX_r)
chrX_m_chrA_r <-subset(heteroz,heteroz[,1] == "chrX_m_chrA_r")
chrX_m_chrX_r <-subset(heteroz,heteroz[,1] == "chrX_m_chrX_r")

# create an empty list that will serve as a container to receive the incoming files
list.data<-list()

# create a loop to read in your data
for (i in 1:length(filelist_bed))
{
  list.data[[i]]<-read.table(filelist_bed[i],sep="\t")
}

# add the names of your data to the list
names(list.data)<-filelist_bed

## Declaration of the dataframes

chr_cases_heteroz <- data.frame()
chr_cases_heteroz <- "NA" # if you don't do that we have error messages
start_position_heteroz <-data.frame()
start_position_heteroz <- "NA"
stop_position_heteroz <-data.frame()
stop_position_heteroz <- "NA"

## Select heteroz sites

# For chrX_m_chrA_r

for (i in 1:length(chrX_m_chrA_r[,1]))
{
  for (j in 1:length(list.data$Xmouse_AND_Arat_ABTC.bed[,2]))
  {
    if(
      (as.numeric(chrX_m_chrA_r[i,2]) > as.numeric(list.data$Xmouse_AND_Arat_ABTC.bed[j,2]) &
       as.numeric(chrX_m_chrA_r[i,2]) < as.numeric(list.data$Xmouse_AND_Arat_ABTC.bed[j,3]) 
       )
      )
    { 
      chr_cases_heteroz[i] <- "chrX_m_chrA_r"
      start_position_heteroz[i] <- as.numeric(list.data$Xmouse_AND_Arat_ABTC.bed[j,2])-1 #we need to substract 1 to have the 1st starting position=0
      stop_position_heteroz[i] <- as.numeric(list.data$Xmouse_AND_Arat_ABTC.bed[j,3])-1
    }
  }
          heteroz <- data.frame(chr_cases_heteroz,
                                start_position_heteroz,
                                stop_position_heteroz)
}

# For chrA_m_chrX_r

for (i in 1:length(chrA_m_chrX_r[,1]))
{
  for (j in 1:length(list.data$Amouse_AND_Xrat_ABTC.bed[,2]))
  {
    if(
      (as.numeric(chrA_m_chrX_r[i,2]) > as.numeric(list.data$Amouse_AND_Xrat_ABTC.bed[j,2]) &
       as.numeric(chrA_m_chrX_r[i,2]) < as.numeric(list.data$Amouse_AND_Xrat_ABTC.bed[j,3]) 
      )
    )
    { 
      chr_cases_heteroz[i] <- "chrA_m_chrX_r"
      start_position_heteroz[i] <- as.numeric(list.data$Amouse_AND_Xrat_ABTC.bed[j,2])-1
      stop_position_heteroz[i] <- as.numeric(list.data$Amouse_AND_Xrat_ABTC.bed[j,3])-1
    }
  }
  heteroz_Am_Xr <- data.frame(chr_cases_heteroz,
                        start_position_heteroz,
                        stop_position_heteroz)
}
# For chrX_m_chrX_r

for (i in 1:length(chrX_m_chrX_r[,1]))
{
  for (j in 1:length(list.data$Xmouse_AND_Xrat_ABTC.bed[,2]))
  {
    if(
      (as.numeric(chrX_m_chrX_r[i,2]) > as.numeric(list.data$Xmouse_AND_Xrat_ABTC.bed[j,2]) &
       as.numeric(chrX_m_chrX_r[i,2]) < as.numeric(list.data$Xmouse_AND_Xrat_ABTC.bed[j,3]) 
      )
    )
    { 
      chr_cases_heteroz[i] <- "chrX_m_chrX_r"
      start_position_heteroz[i] <- as.numeric(list.data$Xmouse_AND_Xrat_ABTC.bed[j,2])-1
      stop_position_heteroz[i] <- as.numeric(list.data$Xmouse_AND_Xrat_ABTC.bed[j,3])-1
    }
  }
  heteroz_Xm_Xr <- data.frame(chr_cases_heteroz,
                        start_position_heteroz,
                        stop_position_heteroz)
}

## Adding everything
ABTC_contigs_start_stop_heteroz_sites <- rbind(unique(heteroz), #need to use 'unique'command of R because we had various heterozygous sites for one contig
                                                      unique(heteroz_Am_Xr),
                                                             unique(heteroz_Xm_Xr))

## Create a document .txt
write.table(ABTC_contigs_start_stop_heteroz_sites,file="chromosomes_ABTC_heteroz.txt",quote=FALSE,col.names=F,row.names=F)


```

OK now I will use vcftools to remove these sections of the RADseq data like this:

``` bash
vcftools --gzvcf MVZ180318-8_concat.fasem_recal_allsites.vcf.gz --exclude-bed ../bedfiles_male_het_sites_to_delete/chromosomes_MVZ_heteroz.bed --out MVZ180318-8_concat.fasem_recal_allsites_minus_malehets.vcf --recode

vcftools --gzvcf ABTC26654-8_concat.fa_recal_allsites.vcf.gz --exclude-bed ../bedfiles_male_het_sites_to_delete/chromosomes_ABTC_heteroz.bed --out ABTC26654-8_concat.fa_recal_allsites_minus_malehets.vcf --recode

```
