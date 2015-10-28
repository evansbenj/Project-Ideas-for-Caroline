# Developing a pipeline for analysis of abyss output

Or go back to the [previous page] (https://github.com/evansbenj/Project-Ideas-for-Caroline/blob/master/Project_9_rat_balls.md).

Caroline and Ben discussed a strategy for making sense of the abyss de novo assembly.  The plan is to take the de novo assembly (filename: ABTC26654-8.fa), which is a fasta file, and blast it against the mouse and rat genomes.  We will save the output as the top hit for each query sequence for each database.  

We think we will start by being 'liberal' with our designation of whether or not a region is on the X chromosome because we can easily test this later by looking for heterozygous sites in males, possibly from the shotgun sequences if the individual is male, which it is not for ABTC26654 but *probably* is for the other two shotgun sequences (check this, and definitely when we map the RADseq data which includes at least one male per species (I think).

We downloaded to the new computer NCBI Blast version 2.2.31+ and we downloaded the mouse and rat genomes we were using (XXX and rn6 respectively) from sharcnet.
