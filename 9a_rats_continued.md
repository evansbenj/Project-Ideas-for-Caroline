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


