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
```
samtools mpileup -d8000 -ugf tropicalis_transcriptome_trinityOut.Trinity.fasta -t DP,AD XT6_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z | bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o XT6_sorted.bam.vcf.gz
```
```
tabix -p vcf XT6_sorted.bam.vcf.gz
```
```
bcftools merge XT2_sorted.bam.vcf.gz XT3_sorted.bam.vcf.gz XT6_sorted.bam.vcf.gz XT9_sorted.bam.vcf.gz XT10_sorted.bam.vcf.gz XT11_sorted.bam.vcf.gz XT16_sorted.bam.vcf.gz XT17_sorted.bam.vcf.gz XT20_sorted.bam.vcf.gz XT1_sorted.bam.vcf.gz XT7_sorted.bam.vcf.gz -Oz XT8_sorted.bam.vcf.gz XT13_sorted.bam.vcf.gz XT19_sorted.bam.vcf.gz -o Merged.vcf.gz
```
```
vcftools --gzvcf Merged.vcf.gz --extract-FORMAT-info AD
```
After renaming, the output file is called Merged.vcf.gz_out.AD.FORMAT and it has exactly what I need - allele depth for all individuals that have genotypes for each transcript for each SNP.

Then I went back to R to get the transcript IDs of significantly female biased and male biased transcripts:

```
borTad_laevisGenome_edgeR_tpm_SL_femalebiased <- borTad_laevisGenome_edgeR_tpm_combine_st46 %>% filter((FDR <= 0.05)& (logFC < -2)) %>% 
  filter(chromosome == "chr8L") %>%
  filter(start<=57000000) 
borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id
write.csv(borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id,
          "borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id.csv", row.names = F)
          
# get the transcript_ids of female_biased trnascripts in the SL region
borTad_laevisGenome_edgeR_tpm_SL_malebiased <- borTad_laevisGenome_edgeR_tpm_combine_st46 %>% filter((FDR <= 0.05)& (logFC > -2)) %>% 
  filter(chromosome == "chr8L") %>%
  filter(start<=57000000) 
borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id
write.csv(borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id,
          "borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id.csv", row.names = F)
```          
And used this list to grep the allele depth data of the individual transcripts from the files generated above with vcftools

```
grep -f SL_trans_IDs_sig_Sex_biased Merged.vcf.gz_out.AD.FORMAT > SL_allelic_depth_sig_sex_biased.AD.FORMAT
```

# Genomic data
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

