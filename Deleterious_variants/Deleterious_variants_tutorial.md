# Deleterious variants analyses

### Introduction

Ensembl's Variant Effect Predictor, __VEP__, is a tool for annotating and analyzing genetic variants, with a focus on deleterious variants. There are prebuild databases for thousands of genomes, and using the tool for these provides some extra features, such as prediction score for deleteriousness of non-synonymous sites, calculated from large sets of homologous proteins. However, it is possible to use the tool on novel genomes as well, using the fasta sequence and a gene annotation file. This will provide one of the following _impact factors_ for all our variants:

| Impact factor | Description |
| --- | ----|
| LOW | synonymous variants |
| MODERATE | non-synonymous variants
| HIGH | non-sense variants (affect start, stop or reading frame) |
| MODIFIER | all other variants (for example intronic, UTRs, intergenic) |

In this tutorial we will first annotate the SNP variants in the Guam rail, and then extract the classes LOW, MODERATE and HIGH to analyze them further.

### Run VEP  

#### 1. Preparing the reference genome and annotation

To run VEP on a non-model organism, you need to have a reference genome and annotation files. The annotation files usually include gene models in GFF or GTF format. The annotation file further needs to be indexed with tabix

Paths to files:
```bash
# assembly file in fasta format
/data/genomics/workshops/smsc_2024/Guam_rail_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta
# Annotation file in gff format
/data/genomics/workshops/smsc_2024/Guam_rail_assembly/Guam_rail.gff
# Variation data for our resequenced individuals in vcf format
/data/genomics/workshops/smsc_2024/???
```
##### a. Create a directory to work in
```bash
mkdir Deleterious
cd Deleterious
```

##### b. Indexing the annotation file
```bash
# Load tabix module
ml ...
# Sort and compress the file
grep -v "#" /data/genomics/workshops/smsc_2024/Guam_rail_assembly/Guam_rail.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c >Guam_rail.gff.gz
# Index
tabix -p gff Guam_rail.gff.gz
```

##### 2. Run the variant predictor
This code should be placed in a script and run on the cluster. If you have your own copy of the assembly file you can use this. You can also make a softlink to the common fasta file using `ln -s originalfile yourfile`.

```bash
# Load ensembl module
ml bio/ensembl-vep
# Run vep
vep -i path/to/variantfile.vcf --cache --gff Guam_rail.gff.gz --fasta bHypOws1_hifiasm.bp.p_ctg.fasta -o vep_rail
```

### Analyze the output

#### 1. Familiarize yourself with the output
If VEP worked, it will produce three output files: vep_rail.txt, vep_rail_summary.html and vep_rail_warnings.txt. The first file contains the raw output. The second file contains a neat summary suitable to open in a web browser (if the hydra cluster doesn't support opening of such file, one option is to save the file locally on your laptop). ADD LINK TO THIS FILE. The third file shows any warnings, for example if there are annotations in the gff file not recognized by VEP.

_For all the commands below, make it a habit to look at the output files, for example with `less`!_

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

#### 2. Extract subsets of data
As we are interested in the deleterious variants, we will mainly focus on the MODERATE and the HIGH impact category. However, to have something potentially neutral to compare with, we will keep the LOW category for a little bit longer.

First we extract the most severe category - nonsense variants
```bash
# Extract variants with a high predicted impact.
grep "IMPACT=HIGH" vep_rail.txt |cut -f2 >high_impact.pos.txt
```
Now do the same for the non-synonymous variants
```bash
grep "IMPACT=MODERATE" vep_rail.txt |cut -f2 >moderate_impact.pos.txt
```
Are there any sites annotated as both HIGH and MODERATE? We can check this by joining the two files:
```bash
join high_impact.variants.txt moderate_impact.pos.txt
```
If there are, we should remove the shared variants from the less severe impact class.
```bash
# Re-run extracting moderate variants, removing positions overlapping with high impact
# join -v1 will return all lines in file 1 not overlapping with file 2.
grep "IMPACT=MODERATE" vep_rail.txt |cut -f2 |join -v1 - high_impact.variants.txt >moderate_impact.pos.txt
```
We can do the same for the LOW impact variants, removing any overlaps with either of the two previous files.
```bash
grep "IMPACT=LOW" vep_rail.txt |cut -f2 |join -v1 - high_impact.pos.txt |join -v1 - moderate_impact.pos.txt >low_impact.pos.txt
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
vcftools -vcf path/to/variantfile.vcf --positions low_impact.pos.txt --max-missing 1 --recode --out low_impact
vcftools -vcf path/to/variantfile.vcf --positions moderate_impact.pos.txt --max-missing 1 --recode --out moderate_impact
vcftools -vcf path/to/variantfile.vcf --positions high_impact.pos.txt --max-missing 1 --recode --out high_impact
```
