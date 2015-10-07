# Background
Any social system can be placed along a axis of sex biased variation in reproductive success. For example, the harem social system of gorillas is thought to have higher variance in male than female reproductive success because a few silverback males dominate most reproduction while other subordinate males do not reproduce whereas most females can reproduce. In contrast, human populations have a more similar variance in reproductive success between the sexes.  The difference in the variance in reproductive success between the sexes has implications for the amount of polymorphism on the sex chromosomes and autosomes. In a system with higher male than female variance in reproductive success, we expect molecular variation on the X chromosome to be elevated and on the Y chromosome to be depressed relative to the expectation with no sex bias in reproductive success.

# Goals of proposed project
Rodents have a striking variation among species in the relative size of the testes to body. Some species have huge testes comprising >3% of the total body mass whereas others have miniscule testes comprising less than 1%. The ratio between testes and body mass is presumably related to the extent of sperm competition within each species. However it is not clear what impact differing levels of sperm competition has on sex biased variance in reproductive success. One possibility is that sperm competition leads to large variance in reproductive success among males.  Alternatively, this could be an example of a "Red Queen" scenario in which competition among males drives the evolution of large testes but without dramatically altering the variance in reproductive success in males.

The goals of the proposed project are to use molecular polymorphism data from rodents with large and small relative testes size to quantify sex biases in the variance in reproductive success with an overall aim of testing how sperm competition does or does not influence reproductive success of males.

# What work is involved
We have two large datasets already in hand including (1) low coverage complete genome sequences for three rodent species that have differing testes:body mass ratios and (2) RADseq polymorphism data for these same species and also a fourth species from multiple individuals.  We have already mapped the RADseq data to the complete genomes of a mouse and rat but failed to identify sufficient X linked markers to adequately test our hypotheses.  As a possible remidy, we hope to assemble low coverage genome sequences by mapping our shotgun sequence reads to the mouse and rat genomes to generate more X-linked regions from our focal species and then map our RADseq data to these regions.  

The proposed steps are as follows:
- for each of the three genomes, map shotgun reads to the mouse and the rat genomes
- examine for each species what proportion of reads map to the X chromosomes of both species; these presumably also are X linked in the focal species
- Do the same for autosomal genes
- Map the RADseq data to putative X linked and autosomal reads and scaffolds for each species
- Test in males whether any putative X linked sites are heterozygous.  If yes, discard these as not X-linked
- Quantify polymorphism and divergence of putative X-linked and autosomal regions in each species; calculate the X:A ratio after correcting for mutation rate

# Anticipated impact
I anticipate this study will be high impact as it will be the first of its kind to test this hypothesis in rodents. We are aware of an analogous study in birds that examines species with sexually dimorphic and non-sexually dimorphic plumage.

# Key collaborators and co-authors
This study was originally conceptialized by Jake Esselstyn and me. Jake is a former postdoc in the lab and now an Assistant Professor at Louisiana State University. Jake and I have co-funded the sequencing and RADseq effort and Jake is a major collaborator on the project.

# Update
Caroline and I are now working on the shotgun genome sequences.  Using `bwa` and `samtools` and a perl script we mapped each of the reads to the mouse genome.  

Here is the perl script:
```perl
#!/usr/bin/perl                                                                                                                                                   
 
use warnings;
use strict;

# This script will execute alignment functions for paired end reads
# launch it in sharcnet like this:
# sqsub -r 7d --mpp=10G -o shotgun.out perl shotgun_alignment.pl

my $path_to_bwa="/work/ben/macaque_RAD_TAGs/bwa-0.6.2";
my $path_to_samtools="/work/ben/macaque_RAD_TAGs/samtools-0.1.18";
my $path_to_data="/work/ben/2015_rat_genomes/ABTC26654";
my $path_to_genome="/work/ben/rat_balls/mouse_genome";
my $genome="mouse_genome_masked.fasta";
my $individual="ABTC26654";

my $commandline;
my $status;

#$commandline = $path_to_bwa."\/bwa aln ".$path_to_genome."\/".$genome." ".$path_to_data."\/".$individual."_R1_trim_paired.fastq \> ".$path_to_data."\/".$individu
al."_R1_trim_paired.fastq.sai";
#$status = system($commandline);
#$commandline = $path_to_bwa."\/bwa aln ".$path_to_genome."\/".$genome." ".$path_to_data."\/".$individual."_R2_trim_paired.fastq \> ".$path_to_data."\/".$individu
al."_R2_trim_paired.fastq.sai";
#$status = system($commandline);
#$commandline = $path_to_bwa."\/bwa aln ".$path_to_genome."\/".$genome." ".$path_to_data."\/".$individual."_R1_trim_single.fastq \> ".$path_to_data."\/".$individu
al."_R1_trim_single.fastq.sai";
#$status = system($commandline);
#$commandline = $path_to_bwa."\/bwa aln ".$path_to_genome."\/".$genome." ".$path_to_data."\/".$individual."_R2_trim_single.fastq \> ".$path_to_data."\/".$individu
al."_R2_trim_single.fastq.sai";
#$status = system($commandline);

$commandline = $path_to_bwa."\/bwa sampe -r \"\@RG\\tID:FLOWCELL1.LANE6\\tSM:".$individual."\\tPL:illumina\" ". $path_to_genome."\/".$genome." ". $path_to_data."\
/".$individual."_R1_trim_paired.fastq.sai ".$path_to_data."\/".$individual."_R2_trim_paired.fastq.sai ".$path_to_data."\/".$individual."_R1_trim_paired.fastq ".$p
ath_to_data."\/".$individual."_R2_trim_paired.fastq \> ".$path_to_data."\/".$individual.".sam";
$status = system($commandline);

$commandline=$path_to_samtools."\/samtools view -bt ".$path_to_genome."\/".$genome." -o ".$path_to_data."\/".$individual.".bam ".$path_to_data."\/".$individual.".
sam";
$status = system($commandline);
$commandline=$path_to_samtools."\/samtools sort ".$path_to_data."\/".$individual.".bam ".$path_to_data."\/".$individual."_sorted";
$status = system($commandline);
$commandline= $path_to_samtools."\/samtools index ".$path_to_data."\/".$individual."_sorted.bam";
$status = system($commandline);
```


Unfortunately only about 2.5% of the reads actually mapped based on the following `samtools` command:
`/work/ben/macaque_RAD_TAGs/samtools-0.1.18/samtools flagstat ABTC26654_sorted.bam`

We thus are trying two alternative strategies:

* de novo assembly using Abyss.  
* direct alignment with stampy


# Abyss de novo assembly

In the `/work/ben/abyss` directory, I have copied a compiled version by someone from sharcnet (Fei), which hopefully will be able to do a de novo assembly.  Before submitting the jobs, you should unload the default modules:

`module unload intel mkl openmpi`

and load the new ones:
```
module load gcc/4.9.2
module load openmpi/gcc-4.9.2/std/1.8.7
module load boost/gcc492-openmpi187/1.59.0
```
and add abyss to the path:
`export PATH=/work/ben/abyss/bin:$PATH`

Then, in this directory: `/work/ben/2015_rat_genomes/ABTC26654`, I entered this command:


`
sqsub -qmpi -n 8 -r 7d --nompirun -o ABTC26654_abyss.out --mpp=10G abyss-pe np=8 name=ABTC26654 k=64 in='ABTC26654_R1_trim_paired.fastq ABTC26654_R2_trim_paired.fastq'
`

If this works, we need to try to optimize the kmer length as described in the Abyss readme file.

# Stampy

* the other approach is to use a program called `stampy` which is specifically designed to map reads to a diverged reference genome. In order to get stampy to work, we need to get the appropriate version of python working on sharcnet. This can be accomplished using the following commands:

```
module unload intel openmpi 
module load gcc/4.8.2 openmpi/gcc/1.8.3 python/gcc/2.7.8 
```

Currently I have formatted the mouse and rat genomes. This was accomplished with the following commands (for mouse):

(first build a genome file)
```
sqsub -r 7d --mpp=4G -o stampy_genome.out ./stampy.py -G mouse_genome_masked /work/ben/rat_balls/mouse_gen
ome/mouse_genome_masked.fasta
sqsub -r 7d --mpp=4G -o stampy_rat_genome.out ./stampy.py -G rn6_masked /work/ben/rat_balls/rat
_genome_rn6/rn6.masked.fa
```
(this will make a genome file with the extension .stidx)

(then build a hash table)
```
sqsub -r 7d --mpp=10G -o stampy_genome_2.out ./stampy.py -g mouse_genome_masked -H mouse_genome_masked
sqsub -r 7d --mpp=10G -o stampy_rat_genome_2.out ./stampy.py -g rn6_masked -H rn6.masked
```
(this will make a file with the extension .sthash)

The mapping of one of the shotgun sequences (ABTC26654) seems to be running.  Here is an example of the sharcnet command I used (within the stampy directory) to do one of the paired end alignment for mouse:

`sqsub -r 7d --mpp=10G -o stampy_genome_3.out ./stampy.py -g ./mouse_genome_masked -h ./mouse_genome_masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_stampy_paired.sam -M /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R1_trim_paired.fastq /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R2_trim_paired.fastq`

Here is an example of the sharcnet command I used (within the stampy directory) to do one of the paired end alignments with rat:

`sqsub -r 7d --mpp=10G -o stampy_genome_3_rat.out ./stampy.py -g ./rn6.masked -h ./rn6.masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_stampy_paired_rat.sam -M /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R1_trim_paired.fastq /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R2_trim_paired.fastq`

## Stampy Update
Unfortunately we hit the 7 day walltime for the trial mouse stampu run.  I am going to start it again on iqaluk using these commands using multiple threads (`-t8`) (MVZ has gz files - much better - cant do the ABTC and JAE uncompressed files yet because abyss is using them):

`/work/ben/stampy-1.0.28/stampy.py -t8 -g ./mouse_genome_masked -h ./mouse_genome_masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_stampy_mouse_paired.sam -M /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_R1_trim_paired.fastq.gz /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_R2_trim_paired.fastq.gz`

and

`/work/ben/stampy-1.0.28/stampy.py -t8 -g ./rn6.masked -h ./rn6.masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_stampy_rat_paired.sam -M /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_R1_trim_paired.fastq.gz /work/ben/2015_rat_genomes/MVZ180318/MVZ180318_R2_trim_paired.fastq.gz`



