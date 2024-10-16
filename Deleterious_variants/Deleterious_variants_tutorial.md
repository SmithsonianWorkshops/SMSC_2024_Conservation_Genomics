# Deleterious variants analyses

## Introduction

Ensembl's Variant Effect Predictor, __[VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)__, is a tool for annotating and analyzing genetic variants, with a focus on deleterious variants. There are prebuild databases for thousands of species, and using the tool on these provides some extra features, such as prediction score for deleteriousness of non-synonymous sites, calculated from large sets of homologous proteins. However, it is possible to use the tool on novel genomes as well, using just the fasta sequence and a gene annotation file. This will provide one of the following _impact factors_ for all our variants:

| Impact factor | Description |
| --- | ----|
| LOW | synonymous variants |
| MODERATE | non-synonymous variants
| HIGH | non-sense variants (affect start, stop or reading frame) |
| MODIFIER | all other variants (for example intronic, UTRs, intergenic) |

We will consider the classes 'MODERATE' and 'HIGH' as potentially deleterious. (This will not always be true, but we will never know the exact effect of all mutations, not even in model organisms, and that's just something we have to live with!)

In this tutorial we will first annotate the SNP variants in the Guam rail, and then extract the classes LOW, MODERATE and HIGH to analyze them further.

## Run VEP  

### 1. Preparing the reference genome and annotation

To run VEP on a non-model organism, you need to have a reference genome and an annotation file. The annotation file usually include gene models in GFF or GTF format. The annotation file further needs to be indexed with tabix.

Paths to files:
```bash
# assembly file in fasta format
/data/genomics/workshops/smsc_2024/Guam_rail_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta
# Annotation file in gff format
/data/genomics/workshops/smsc_2024/Guam_rail_assembly/Guam_rail.gff
# Variation data for our resequenced individuals in vcf format
/data/genomics/workshops/smsc_2024/VCF/gatk_GuamRails_autosomes_header.vcf
```
#### a) Create a directory to work in
```bash
mkdir Deleterious
cd Deleterious
```

#### b) Indexing the annotation file
Vep requires that exons are annotated. Since our annotation file lacks exons, we will duplicate the CDS annotation and just replace the "CDS" with "exon" before we index the file:
```bash
# Tabix should be included in the vep module
module load bio/ensembl-vep
# Sort file and add exons before compressing
grep -v "#" /data/genomics/workshops/smsc_2024/Guam_rail_assembly/Guam_rail.gff | \
sort -k1,1 -k4,4n -k5,5n -t$'\t' | \
awk -v OFS="\t" '{if($3=="CDS"){line=$0; $3="exon"; $8="."; print line"\n"$0}else{print $0}}' | \
bgzip -c >Guam_rail.withExons.gff.gz
# Index
tabix -p gff Guam_rail.withExons.gff.gz
```

#### c) Prepare the vcf file
If we want we could just use the vcf file as it is, but we can also do some more filtering. Here I decided to remove indels, and only keep bi-allelic sites (SNPs with two alleles, not more) with no missing data:
```bash
ml bio/vcftools/0.1.16
vcftools --vcf gatk_GuamRails_autosomes_header.vcf --remove-indels  \
 --max-missing 1.0 --max-alleles 2 --recode --out GuamRails_SNPs
```

### 2. Run the variant predictor
This code should be placed in a script and run as a job. It took around ~30 min using 8 cores.  If you already have your own copy of the assembly file you can use that. You can also make a softlink to the common fasta file using `ln -s originalfile yourfile`.

```bash
# Load ensembl module
module load bio/ensembl-vep
# Run vep
vep -i GuamRails_SNPs.recode.vcf --fork 8 --offline --gff Guam_rail.withExons.gff.gz \
 --fasta bHypOws1_hifiasm.bp.p_ctg.fasta -o vep_rail.txt --species gallus_gallus --force_overwrite
#--fork 8 tells vep to use 8 cores - make sure to give it as many cores when you start the script.
#--offline means we don't want to use the vep database but our own data
#--species gallus_gallus means... well, we know it's not a chicken, but without this flag
# vep assumes a human genome, and giving it a custom name returned an error when I tested the code.
```

## Analyze the output

### 1. Familiarize yourself with the output
If VEP worked, it will produce two or three output files: vep_rail.txt, vep_rail.txt_summary.html and possibly vep_rail.txt_warnings.txt. The first file contains the raw output. The second file contains a neat summary suitable to open in a web browser (if the hydra cluster doesn't support opening this file, one option is to save the file locally on your laptop). (A copy of the summary file can be found [HERE](https://github.com/SmithsonianWorkshops/SMSC_2024_Conservation_Genomics/blob/main/Deleterious_variants/vep_GuamRails_summary.html)). The third file shows any warnings, for example if there are annotations in the gff file not recognized by VEP.

_For all the commands below, make it a habit of looking at the output files, for example with `less`!_

All the information in the summary can also be found in the raw output. We can use a variety of unix tools to check it out:
```bash
# Count the number of annotated variants
# (Grep -v removes the header lines: the ones that start with #)
grep -v "^#" vep_rail.txt |wc -l
```
How does this relate to the number variants in your input file?
```bash
# Count the number of unique annotated variants
# (cut -f1 extracts the first column only, uniq take only unique lines, and wc -l counts the lines)
grep -v "^#" vep_rail.txt |cut -f1 |uniq |wc -l
```

Many variants are annotated several times, for example if a gene has multiple transcripts, or if two genes are overlapping.

```bash
# Check what annotations are classified as having a low impact, and how many there are of each type
# (sort will sort the output in alphabetical order, and uniq -c will count all unique lines)
grep "IMPACT=LOW" vep_rail.txt  |cut -f7 |sort |uniq -c
```

Now you can do the same with for the other impact types.

### 2. Extract subsets of data
As we are interested in the deleterious variants, we will mainly focus on the MODERATE and the HIGH impact category. However, to have something potentially neutral to compare with, we will keep the LOW category for a little bit longer.

First we extract the most severe category - nonsense variants
```bash
# Extract variants with a high predicted impact.
grep "IMPACT=HIGH" vep_rail.txt |cut -f2 |uniq >high_impact.pos.txt
```
Now do the same for the non-synonymous variants
```bash
grep "IMPACT=MODERATE" vep_rail.txt |cut -f2 |uniq >moderate_impact.pos.txt
```
Are there any sites annotated as both HIGH and MODERATE? We can check this by joining the two files:
```bash
join high_impact.variants.txt moderate_impact.pos.txt
```
If there are, we should remove the shared variants from the less severe impact class.
```bash
# Re-run extracting moderate variants, removing positions overlapping with high impact
# join -v1 will return all lines in file 1 not overlapping with file 2.
grep "IMPACT=MODERATE" vep_rail.txt |cut -f2 |uniq |join -v1 - high_impact.variants.txt >moderate_impact.pos.txt
```
We can do the same for the LOW impact variants, removing any overlaps with either of the two previous files.
```bash
grep "IMPACT=LOW" vep_rail.txt |cut -f2 |uniq |join -v1 - high_impact.pos.txt |join -v1 - moderate_impact.pos.txt >low_impact.pos.txt
```

If we want to extract the variants from the vcf file using vcftools, we need tab separated positions files. We can use the Stream EDitor sed to replace the ":" to tabs (\t).
```bash
sed -i '' 's/:/\t/g' low_impact.pos.txt
sed -i '' 's/:/\t/g' moderate_impact.pos.txt
sed -i '' 's/:/\t/g' high_impact.pos.txt
```

Create new vcf files with the variants we are interested in.
```bash
ml bio/vcftools/0.1.16
vcftools --vcf GuamRails_SNPs.recode.vcf --positions low_impact.pos.txt \
--recode --out GuamRail_low_impact
vcftools --vcf GuamRails_SNPs.recode.vcf --positions moderate_impact.pos.txt \
--recode --out GuamRail_moderate_impact
vcftools --vcf GuamRails_SNPs.recode.vcf --positions high_impact.pos.txt \
--recode --out GuamRail_high_impact
```

### 3. Analyze the alleles
Before we start investigating genetic load, there is one important thing we need to think about: _Which_ of the two alleles in a site is the deleterious one?? Vep doesn't provide us with this information. We know that a mutation in a certain position causes a non-synonymous variation, and that this could be harmful. But for all we know, it could be our reference individual who is carrying the deleterious allele, and the other individuals carrying the 'normal' variant.

There are a few different ways to figure this out (see further below*), but for now we will assume that the REFERENCE allele is 'normal', and that the ALTERNATIVE allele is deleterious.

A very convenient tool to count allele types and plot results is the vcfR package in R. As this course is not focusing on R, I'll show how we can look at genotypes and count alleles using different unix tools and vcftools. For those interested in the vcfR code, it can be found at the bottom of this page**.

#### a) Allele frequency spectrum
A good way to see if our potentially deleterious sites are under more selective constraints than for example synonymous mutations, is to compare their allele frequency spectra. With only two individuals we only have four alleles to work with, but let's give it a try!

```bash
# First we'll just look at the genotypes (remove everything else)
grep -v "##" GuamRail_moderate_impact.recode.vcf | \
 awk '{out=$1"\t"$2; for(i=10; i<=11; i++){split($i,s,":"); out=out"\t"s[1]}; print out}' |less
```
This is a good start to just get a feeling for the data.

Now we will use vcftools to calculate allele frequency for each site. This will create .frq output files that we can summarize into a little table.
``` bash
# Create a table header
echo "Type Number Ref_freq Alt_freq" |sed 's/ /\t/g' >SFS.txt
# Loop over the types, create a frequency table and summarize
for type in "low" "moderate" "high"
do
  vcftools --vcf GuamRail_${type}_impact.recode.vcf --freq --out GuamRail_${type}_impact
  tail -n+2 GuamRail_${type}_impact.frq  |cut -f5,6 |sed 's/:/\t/g' | \
   cut -f2,4 |sort |uniq -c |awk -v t=$type -v OFS="\t" '{print t,$0}' >>SFS.txt
done
# A lot of different unix tools here, including an awk script..
# Make sure you know what each step does!
```
The file SFS.txt contains the site frequency spectra for all the three types of mutations. Below is some R code you can use for plotting, but you can use any tool you like (even Excel)

```R
require(tidyverse)

file<-"SFS.txt"
SFS<-file %>% read.table(header=TRUE) %>% as_tibble()

# Order the types
SFS$Type<-factor(SFS$Type, levels=c("high","moderate","low"))

ggplot(SFS, aes(x=Alt_freq, y=Number, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```
What do we see here? The sizes are so different between the types so they are hard to compare! We can try making the bars using relative sizes instead:

```R
# With relative sizes
SFS_rel <- SFS %>% group_by(Type) %>% mutate(Rel=Number/sum(Number))

ggplot(SFS_rel, aes(x=Alt_freq, y=Rel, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```
Now it looks better! We see that the 'HIGH' category is shifted to the left, and the 'LOW' category has relatively more sites with higher alternative frequency. Can you explain why?


#### b) Masked and realized load
From the lectures, we recall that masked load comes from heterozygous deleterious mutations, and realized load comes from homozygous deleterious mutations. First we can just look at different genotype counts from one individual:

```bash
# First individual in column 10 (cut -f10 means extract the 10th column)
grep -v "##" GuamRail_moderate_impact.recode.vcf |cut -f10 |cut -f1 -d":" |sort |uniq -c
```
To look at the next individual we can replace `cut -f10` with `cut -f11`. If we haven't noticed it before, some sites are _phased_, with a '|' between the alleles instead of the standard '/'. When we count for example heterozygous sites, we take both "0/1" and "0|1". Do you notice some striking difference between the individuals?

Now we use some more awk and unix tools to save heterozygous and homozygous alternative sites from high and moderate separately
```bash
# Create a new file with just a header
echo "Type Ind Load Number" |sed 's/ /\t/g' >Load_table.txt
# Loop over the types, create a frequency table and summarize
for type in "moderate" "high"
do
  # Loop over the individuals (in column 10 and 11)
  for col in {10..11}
  do
    # Save the name of the individual
    ind=`grep "#CHR" GuamRail_${type}_impact.recode.vcf |cut -f $col`
    # Extract the genotypes and use awk to count the heterozygous and
    # homozygous as masked and realized respctively
    grep -v "#" GuamRail_${type}_impact.recode.vcf |cut -f$col |cut -f1 -d":" | \
     awk -v OFS="\t" -v t=$type -v i=$ind -v masked=0 -v realized=0 '{if($1=="0/1" || $1=="0|1"){masked++}else if($1=="1/1" || $1=="1|1"){realized++}}END{print t,i,"masked",masked,"\n"t,i,"realized",realized}' >>Load_table.txt
  done
done
```

Again, we can plot this with R:
```R
require(tidyverse)
file<-"Load_table.txt"
LOAD<-file %>% read.table(header=TRUE) %>% as_tibble()

# Make one bar for each individual, and one plot for each mutation type
ggplot(LOAD, aes(x=Load, y=Number, fill=Ind)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.8) +
  facet_wrap(Type~., nrow=2, scales='free_y')
```
Do you see a difference between the individuals? Or a difference between 'HIGH' impact mutations and 'MODERATE' impact mutations?

Remenber, this is just a very rough estimation of genetic load! Apart from finding the correct deleterious allele (which we ignored above), it might be necessary to account for differences in sequencing (if some individuals have more missing data, for example). One way to do this is to calculate load as the number of deleterious alleles per genotyped sites.  


This is the end of this tutorial! For the interested there are some extra information and code below.

---------
### *About finding out which allele is deleterious
In a large population, deleterious variants will most likely be selected against, and will never reach high frequencies. Therefore it is often safe to assume that the _minor_ allele is the deleterious variant. But in a small population, we know that also deleterious variants can reach high frequencies just due to drift. And what if we only have a couple of samples in the first place! With only 4 alleles in total, it is hard to tell which is the minor!
Another option is to _polarize_ the variants (i.e. to determine the ancestral allele), and assume the ancestral allele is the 'normal' variant. This approach has been used in many conservation studies. Can you think of caveats with this method?

Do we have enough data in the Guam Rail project to polarize our variants?

---------
### ** Repeat Analysis part 3. in R with the vcfR package
This code is adapted from the github [https://github.com/linneas/wolf-deleterious](https://github.com/linneas/wolf-deleterious) (Smeds et al. 2023).
When using vcfR it is more convenient to have all variants in a single vcf file:
```bash
rm -f all_three.pos.txt
for type in "low" "moderate" "high"
do
  awk -v t=$type '{print $0"\t"t}' ${type}_impact.pos.txt |uniq >>all_three.pos.txt
done
vcftools --vcf GuamRails_SNPs.recode.vcf --positions all_three.pos.txt \
--recode --out GuamRail_all_three
```
Read the vcf into R and perform the same analysis as above:
```R

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)

# Reading the two files
vcf_file="GuamRail_all_three.recode.vcf"
type_file="all_three.pos.txt"
vcf <- read.vcfR(vcf_file)
vep_types <- type_file %>% read.table(header=FALSE) %>% as_tibble() %>%
        rename(CHROM=V1, POS=V2, Type=V3)

# Convert vcf to tidy format
tidy_vcf <- vcfR2tidy(vcf,
              #info_fields=c("AA"),   # This could be used if we have ancestral information
              format_fields=c("GT"),
              dot_is_NA=TRUE)

# Extracting relevant data from vcf and merge with VEP impact classes
tidy_gt <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% inner_join(tidy_vcf$fix) %>%
  inner_join(vep_types) %>% mutate(gt=paste(substr(gt_GT,1,1),substr(gt_GT,3,3), sep=""))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Building the SFS
SFS <- tidy_gt %>% select(CHROM, POS, Indiv, gt, Type) %>%
  group_by(CHROM, POS, gt, Type) %>% summarize(count=n()) %>%
  mutate(alt=case_when((gt=='11') ~ (count*2.0),
                       (gt=='01') ~ (count*1.0),
                       TRUE ~ 0)) %>%
  group_by(CHROM, POS, Type) %>%  summarize(totalt=sum(alt)) %>%
  group_by(Type, totalt) %>% summarize(count=n()) %>% ungroup() %>%
  group_by(Type) %>% mutate(Rel=count/sum(count))

# Order the types before plotting
SFS$Type<-factor(SFS$Type, levels=c("high","moderate","low"))

# Plot SFS
ggplot(SFS, aes(x=totalt, y=Rel, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(x="Alt alleles", y="Fraction of sites")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate the load

# Summarize the number of genotypes (skip 'low' mutations)
LOAD <- tidy_gt %>% filter(Type!="low") %>% select(Indiv, Type, gt) %>%
            group_by(Indiv,Type,gt) %>% summarize(count=n()) %>% ungroup() %>%
            mutate(Load=case_when(gt=="01" ~ "masked",
                                  gt=="11" ~ "realized",
                                  TRUE ~ NA)) %>% drop_na()

# Plot load
ggplot(LOAD, aes(x=Load, y=count, fill=Indiv)) +
  geom_bar(stat="identity", position=position_dodge(), alpha=0.8) +
  facet_wrap(Type~., nrow=2, scales='free_y')
```
