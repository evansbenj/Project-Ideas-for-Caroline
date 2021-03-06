# Background
A recent [study](http://www.pnas.org/content/112/34/E4752) based on laboratory animals has reported polymorphism in the sex chromosomes of the frog *X. tropicalis*.  Because this study is based on lab animals of unknown provenance, it is not clear whether this variation in sex chromosomes naturally occurs within one population or whether it occurs in a geographically structured fashion

# Goals
To evaluate whether there is geographic structure to the mechanism of sex determination in the frog *X. tropicalis*

# Work to be done
This project has not yet been started.  It will involve fieldwork in Ghana which is scheduled for mid-year 2016.  On this expedition, we hopefully can collect *X. tropicalis* from multiple localities and export live aniamls as well as preserved animals. These animals could then be crossed in the lab, raised to postmetamorphosis, and screened for ZW versus XY sex systems using GBS.

# Anticipated Impact
I anticipate that this work will result in a medium to high impact publication suitable for a journal such as Evolution or Genome Biology and Evolution.  *Xenopus tropicalis* is a widely used organism for research and the genetics behind sex determination would thus be of keen interest to many.  This project would expose Caroline to fieldwork in Africa, and GBS bioinformatic approaches using a reference genome.

# Key collaborators
This would be a collaborative effort in the Evans lab with possible outside contributors, possibly led by Caroline.  This project cannot begin until fieldwork is completed in mid 2016, so other projects will be started in the meantime.

# Update from June, 2020
With Caroline and BenF, we did fieldwork in Ghana and brought back frogs which turned out to have Y chromosomes.  Now we are revising a paper and need to assess evidence for degeneration of the W-chromosome and also the potential for reference bias in the detection of male biased transcripts.  We are using three approaches for this:
* If Y is derived from Z and females have degeneration of the W, then males would be expected to have higher expression. 
* If reference bias is an issue - level of female expression should be lower than the autosomal. We can compare the median expression of SL and non-SL transcripts in females. Maybe get the ratio of these expression levels by individual.
* Map RNAseq data to v10 or the transcriptome assembly. Count SNP/bps of non-male-specific male-biased transcripts in SL region and compare to SNPs/bp of non-male biased transcripts in SL region.

# RNAseq
RNAseq data for trop tads is here:
```
/home/evanslab/trop_tadpole_RNAseq/data/trimmed_RNAseq_data
```
Transcriptome assembly is here:
```
/home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta
```
Here are some of the files that Xue generated with a subset of individuals:
Group1 and group2 corresponding to analysis 1 and 2 respectively in the previous email.

Group 1 has these WW females: XT3, XT9, XT11, and XT20 and these WY males: XT7 and XT8

Group 2 has these WZ females: XT2, XT6, XT10, XT16, and XT17 and these ZY males: XT1, XT13, XT19

Group 2 has these WW females: XT3, XT9, XT11, and XT20 and these ZY males: XT1, XT13, XT19
```
/home/evanslab/trop_tadpole_RNAseq/data/combine_count_de_detail/expression_tpm/trop_tad_v10_tpm_uncollapsed_group1_combine.csv
/home/evanslab/trop_tadpole_RNAseq/data/combine_count_de_detail/expression_tpm/trop_tad_v10_tpm_uncollapsed_group2_combine.csv
/home/evanslab/trop_tadpole_RNAseq/data/combine_count_de_detail/expression_tpm/trop_tad_v10_tpm_uncollapsed_group3_combine.csv
```
I'm mapping the RNAseq data to the transcriptomes here:
```
/home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome
```
Using these commands:
```
bwa mem tropicalis_transcriptome_trinityOut.Trinity.fasta '<zcat /home/evanslab/trop_tadpole_RNAseq/data/trimmed_RNAseq_data/XT6_R1_paired.fastq.gz' '<zcat /home/evanslab/trop_tadpole_RNAseq/data/trimmed_RNAseq_data/XT6_R2_paired.fastq.gz' > XT6.sam
```
```
samtools view -bt tropicalis_transcriptome_trinityOut.Trinity.fasta -o XT6.bam XT6.sam
```
```
rm -f XT6.sam
```
```
samtools sort XT6.bam -o XT6_sorted.bam
```
```
samtools index XT6_sorted.bam
```
Call genotypes, do light filter, output all calls, including homoz REF
```
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT10_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT10_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT11_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT11_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT13_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT13_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT16_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT16_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT17_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT17_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT19_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT19_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT1_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT1_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT20_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT20_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT2_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT2_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT3_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT3_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT6_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT6_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT7_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT7_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT8_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT8_sorted.bam.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT9_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XT9_sorted.bam.vcf.gz
```
```
tabix -p vcf XT6_sorted.bam.vcf.gz
```
```
bcftools merge XT2_sorted.bam.vcf.gz XT3_sorted.bam.vcf.gz XT6_sorted.bam.vcf.gz XT9_sorted.bam.vcf.gz XT10_sorted.bam.vcf.gz XT11_sorted.bam.vcf.gz XT16_sorted.bam.vcf.gz XT17_sorted.bam.vcf.gz XT20_sorted.bam.vcf.gz XT1_sorted.bam.vcf.gz XT7_sorted.bam.vcf.gz -Oz XT8_sorted.bam.vcf.gz XT13_sorted.bam.vcf.gz XT19_sorted.bam.vcf.gz -o Merged.vcf.gz
```

I made this into a tab file and I'm running my polymorphism script on it. There are three possible sex chromosome genotypes for our tad family:
* WZ mother x WY father
* WW mother x WY father
* WW mother * ZY father.

For the first option, daughters are WW or WZ and sons are WY or ZY.  On average, if the XT Y is from a Z and the W is degenerate, I expect more polymorphism in sex-linked transcripts from males than females in this scenario, especially in the male-biased transcripts. Polymorphism in male-biased sex-linked transcripts should be lower than polymorphism in male-biased non-sex-linked transcripts.

For the second option, if the XT Y is from a Z and the W is degenerate, I expect more polymorphism in sex-linked transcripts from males than females in this scenario. Polymorphism in male-biased sex-linked transcripts should be lower than but fairly similar to polymorphism in male-biased non-sex-linked transcripts.

For the third option, if the XT Y is from a Z and the W is degenerate, I'd expect a similar level of polymorphism in the males and females. I would also not expect such a male-biased expression skew, so this cross a priori is unlikely.

Probably can't distinguish between the first and second options.

However, if the Y were derived from a degenrate W, the predictions would be different.

The first option: similar polymorphism expected in female sex-linked transcripts (WW or WZ) as male sex-linked transcripts (WY or ZY)
The second option: similar polymorphism expected in female sex-linked transcripts (WW) as male sex-linked transcripts (WY)
The third option: more polymorphism expected in female sex-linked transcripts (WZ) than male sex-linked transcripts (WY)

```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts_.pl XT_st50_RNAseq_Merged.vcf.gz.tab 11111111111111 3_4_1_2_3_4_5_6_7_8_9_10_11_12_13_14 XT_tads_RNAseq_Merged_females_and_males_.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts_.pl XT_st50_RNAseq_Merged.vcf.gz.tab 11111111111111 3_4_1_2_4_5_8_9_10_11_14 XT_tads_RNAseq_Merged_females_.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts_.pl XT_st50_RNAseq_Merged.vcf.gz.tab 11111111111111 3_4_3_6_7_12_13 XT_tads_RNAseq_Merged_males_.poly_by_windows
```

```
vcftools --gzvcf Merged.vcf.gz --extract-FORMAT-info AD
```
After renaming, the output file is called Merged.vcf.gz_out.AD.FORMAT and it has exactly what I need - allele depth for all individuals that have genotypes for each transcript for each SNP.

I thought that I identified a bug that caused the AD numbers to be weird from vcftools and bcftools (but it was excel; see below). So I wrote my own script to get this:
```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file and pulls out allele depth information.  This is possible with vcftools
# and bcftools but those programs both make errors.

# to execute type 
# ./gets_AD_from_Vcf_filter.pl Merged.vcf.gz 14 Merged.vcf.gz.out.AD.FORMAT
# where number_of_samples refers to the number_of_samples in the vcf file

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[2];


if ($inputfile =~ /.gz$/) {
open(DATAINPUT, "gunzip -c $inputfile |") || die "can’t open pipe to $inputfile";
}
else {
open(DATAINPUT, $inputfile) || die "can’t open pipe to $inputfile";
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my $number_of_samples = $ARGV[1];

my @male_chrX_hets=();
my @male_chrY_hets=();
my @male_chrX_lowcoverage=();
my @female_chrX_lowcoverage=();
my @male_chrY_lowcoverage=();
my @female_chrY=();
my @autosomal_lowcoverage=();
my @male_Yhet_sites;
my @male_Xhet_sites;
my @female_Y_sites;
my $y;
my @columns=();
my @DP;
my $AD;
my @genotypes;

while ( my $line = <DATAINPUT>) {
	# print all commented lines to the outfile
	if(substr($line,0,2) eq "#C"){
		print OUTFILE $line;
	}
	elsif(substr($line,0,1) ne "#"){
		@columns=split(/\s/,$line);
		# first find out where the depth of coverage is			
		@DP = split(":",$columns[8]);
		$AD=-1;
		for ($y=0; $y<= $#DP; $y += 1){
			if($DP[$y] eq "AD"){
				$AD=$y;
			}
		}
		# now print the first 5 columns
		for ($y=0; $y<= 4; $y ++){
			print OUTFILE $columns[$y],"\t";
		}
		# now print out the Allele Depth of each sample
		for ($y=9; $y<= ($number_of_samples + 8); $y ++){
			# check each individual genotype in the ingroup
			@genotypes=split(":",$columns[$y]);
			# if a call was made, check if either genotupe is more than one base
			print OUTFILE $genotypes[$AD],"\t";
		} # end for loop over each genotype within a base position
		print OUTFILE "\n";
	} # end else	
}# end while
close DATAINPUT;
close OUTFILE;
````
Turns out it was actually an issue with the way Excel imported the text file! Anyhow, doesn't matter.

# Now get transcript IDs of male-biased and female-biased trnascripts in SL region
```R
XT_data_whole_genome_tpm_SL_malebiased <- XT_data_whole_genome %>% 
  filter((FDR <= 0.05)& (logFC > 2)) %>% 
  filter(chromosome == "Chr7") %>%
  filter(start<=12000000) 
XT_data_whole_genome_tpm_SL_malebiased$trans_id
write.csv(XT_data_whole_genome_tpm_SL_malebiased$trans_id,
          "XT_data_whole_genome_tpm_SL_malebiased$trans_id.csv", row.names = F)

XT_data_whole_genome_tpm_SL_femalebiased <- XT_data_whole_genome %>% 
  filter((FDR <= 0.05)& (logFC < -2)) %>% 
  filter(chromosome == "Chr7") %>%
  filter(start<=12000000) 
XT_data_whole_genome_tpm_SL_femalebiased$trans_id
write.csv(XT_data_whole_genome_tpm_SL_femalebiased$trans_id,
          "XT_data_whole_genome_tpm_SL_femalebiased$trans_id.csv", row.names = F)
```
And use this to grep the allele depth data of the individual transcripts from the file generated above with vcftools
```bash
grep -f female_biased_SL_transIDs.txt Merged.vcf.gz.out.AD.FORMAT > SL_allelic_depth_sig_female_biased.AD.FORMAT.txt
grep -f male_biased_SL_transIDs.txt Merged.vcf.gz.out.AD.FORMAT > SL_allelic_depth_sig_male_biased.AD.FORMAT.txt
```
this above stuff was done here:
```
/home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome
```

# RADseq data
I'm working with some old RADseq vcf files mapped to v9 by BenF here: `/home/evanslab/tropicalis_families`

```
vcftools --gzvcf mpileup_raw_tropicalisFamilyEast_allChrs.vcf.gz --indv BJE4687_girl_GhE.fq.gz --chr Chr07 --window-pi 1000000
```
For Ghana East, males and females both have peaks <20Mb on chr7.  THis is consistent with a WWxYZ cross, with daughters being WZ and sons being WY.

# Genomic data

Using bam file from Berkeley on graham: `/home/ben/scratch/XT-v10`
```
samtools depth -r Chr7:1-12000000 JBL052.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)*(sum/NR))}'
Average =  54.3531
Stdev =  65.6165
```
```
samtools depth -r Chr7:12000001-133565930 JBL052.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)*(sum/NR))}'
Average =  60.0115
Stdev =  53.3348
```

Looked at diversity:
```
vcftools --gzvcf JBL052_chr7.bam.vcf.gz --window-pi 100000
```
It is much higher on Chr7 <16.5 Mb than the rest of the chr.  This is a ZW individual.

Trying to figure out depth of XT reads in v10.  Only have shitty 454 seqs which I mapped here:
```
/4/ben/XT_v10
```
Coverage is ridiculously low but I can check depth on the covered bases and it is about the same in and out of the SL region:

```
[ben@info115 XT_v10]$ samtools depth  SRR000914_to_v10_sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
Average =  1.17029
[ben@info115 XT_v10]$ samtools depth -r Chr7:1-12000000 SRR000914_to_v10_sorted.bam  |  wc -l
596535
[ben@info115 XT_v10]$ samtools depth  SRR000914_to_v10_sorted.bam  | wc -l
100776206
[ben@info115 XT_v10]$ pwd
```

Another idea is to map assembled transcripts to v10 and quantify how many divergent sites are in sex-linked, male-biased and other transcripts. If the Y came from the Z, we'd expect more divergent sites in SL than nonSL transcripts if the ref is a WW individual.  And this may especially be the case for SL male biased transcipts.

I installed a tool to tabulate divergence from ref using a bam file here (http://lindenb.github.io/jvarkit/Sam2Tsv.html):
```
/home/evanslab/trop_tadpole_RNAseq/data/jvarkit/dist
```
I'm working with a bam file that Xue made using gmap, which is splice aware, here:
```
/home/evanslab/trop_tadpole_RNAseq/data/mapping_tropDNTrans_tropGenomeV10_gmap
```
After using picard CreatDictionary and samtools faidx on the ref, I made a tab delimited file like this:
```
java -jar /home/evanslab/trop_tadpole_RNAseq/data/jvarkit/dist/sam2tsv.jar -r /home/evanslab/trop_tadpole_RNAseq/data/XT_v10/XENTR_10.0_genome.fasta.gz /home/evanslab/trop_tadpole_RNAseq/data/mapping_tropDNTrans_tropGenomeV10_gmap/tropDNTtrans_tropGenomeV10_gmap.bam | perl -ane 'print if ($F[4] ne uc($F[7]))&&($F[4] ne ".")&&($F[4] ne "*")' > divergent_sites.txt
```
This pipes the output to perl to parse and save only the divergent sites that are not gaps in the transcripts or asterisks, which may mean that they were hard clipped.

The output filee is huge, so I am going to filter it to focus only on chr7, and include all matching sites, not only the mismatching ones. This is important because we need to know the propotion of the alignment that is diverged.
```
java -jar /home/evanslab/trop_tadpole_RNAseq/data/jvarkit/dist/sam2tsv.jar -A -r /home/evanslab/trop_tadpole_RNAseq/data/XT_v10/XENTR_10.0_genome.fasta.gz /home/evanslab/trop_tadpole_RNAseq/data/mapping_tropDNTrans_tropGenomeV10_gmap/tropDNTtrans_tropGenomeV10_gmap.bam | perl -ane 'print if ($F[3] eq "Chr7")&&($F[4] ne ".")&&($F[4] ne "*")&&($F[8] eq "M")' > divergent_sites_A.txt
```
I took out only Chr7 and then ran the script like this:
```
samtools sort tropDNTtrans_tropGenomeV10_gmap.bam -o tropDNTtrans_tropGenomeV10_gmap_sorted.bam
```
```
samtools index tropDNTtrans_tropGenomeV10_gmap_sorted.bam
```
```
samtools view -b tropDNTtrans_tropGenomeV10_gmap_sorted.bam Chr7 > tropDNTtrans_tropGenomeV10_gmap_sorted_Chr7only.bam
```
```
java -jar /home/evanslab/trop_tadpole_RNAseq/data/jvarkit/dist/sam2tsv.jar -A -r /home/evanslab/trop_tadpole_RNAseq/data/XT_v10/XENTR_10.0_genome.fasta.gz /home/evanslab/trop_tadpole_RNAseq/data/mapping_tropDNTrans_tropGenomeV10_gmap/tropDNTtrans_tropGenomeV10_gmap_sorted_Chr7only.bam | perl -ane 'print if ($F[4] ne ".")&&($F[4] ne "*")&&($F[8] eq "M")' > divergent_sites_A.txt
```

# Extract sequences from transcriptome assembly

Candidates for kos (/home/evanslab/trop_tadpole_RNAseq/data/build_transcriptome):

```
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(TRINITY_DN1430_c1_g2_i3 TRINITY_DN20482_c0_g1_i4 TRINITY_DN5596_c0_g1_i3 TRINITY_DN21880_c0_g1_i1 TRINITY_DN12616_c0_g1_i1 TRINITY_DN23111_c0_g1_i1 TRINITY_DN2118_c0_g1_i8 TRINITY_DN21513_c0_g1_i6)}print if $c' tropicalis_transcriptome_trinityOut.Trinity.fasta  >> KO_candidates.fa
```
