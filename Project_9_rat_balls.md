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
Caroline and I are now working on the shotgun genome sequences.  Using `bwa` and `samtools` and a perl script we mapped each of the reads to the mouse genome.  Unfortunately only about 2.5% of the reads actually mapped based on the `samtools faidx` command.  We thus are trying to alternative strategies:
* de novo assembly using Abyss.  I "may" have managed to install this as follows:
```
autoreconf
./configure --prefix=/work/ben/abyss-1.9.0
automake
make
```
The autoreconf and automake stuff was needed to avoid an error `aclocal-1.13: command not found`.

In the `/work/ben/abyss-1.9.0/bin` directory, I have entered the following command, which hopefully will start a de novo assembly: 
`sqsub -r 7d --mpp=10G -o abyss_ABTC26654_mouse.out ./abyss-pe name=ABTC26654_mouse k=64 in='/work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R1_trim_paired.fastq /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R2_trim_paired.fastq'`

If this works, we need to try to optimize the kmer length as described in the Abyss readme file.

* the other approach is to use a program called `stampy` which is specifically designed to map reads to a diverged reference genome. Currently I have formatted the mouse and rat genomes and the mapping of one of the shotgun sequences (ABTC26654) seems to be running.  Here is an example of the sharcnet command I used (within the stampy directory) to do one of the paired end alignment for mouse:
`sqsub -r 7d --mpp=10G -o stampy_genome_3.out ./stampy.py -g ./mouse_genome_masked -h ./mouse_genome_masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_stampy_paired.sam -M /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R1_trim_paired.fastq /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R2_trim_paired.fastq`

Here is an example of the sharcnet command I used (within the stampy directory) to do one of the paired end alignments with rat:
`sqsub -r 7d --mpp=10G -o stampy_genome_3_rat.out ./stampy.py -g ./rn6.masked -h ./rn6.masked  --substitutionrate=0.10 -o /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_stampy_paired_rat.sam -M /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R1_trim_paired.fastq /work/ben/2015_rat_genomes/ABTC26654/ABTC26654_R2_trim_paired.fastq`

 

