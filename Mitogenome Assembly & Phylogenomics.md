# Mitogenome Assembly & Phylogenomic Tutorial

Heterozygosity is a simple and informative statistic that can be obtained by analyzing whole-genome data. You can calculate the average heterozygosity of an individual or assess the local heterozygosity throughout the genome to answer more specific questions.

There are many methods available to achieve this. In this tutorial, we will use two strategies. To obtain the average heterozygosity, we will manipulate the VCF file using bcftools and vcftools. To explore heterozygosity variation across the genome, we will use ANGSD and other associated tools.

### Mitogenome Assembly using MitoZ

First we want to subsample the whole genome reads for each individual we have. For the purposes of practicing we are only going to use our two resequenced individuals for the mitogenome assembly but for the whole phylogeny later we are going to have different individuals included

```
mkdir mitogenomes/
mkdir mitogenomes/sub_reads
cd /mitogenomes/
```

```bash
#!/bin/bash
# randomly subsample 20% of read pairs
seqtk sample -s 100 /data/genomics/workshops/smsc_2024/rawdata/HOW_N23-0063_1.fastq.gz  0.2 | pigz -p 8 > sub_reads/HOW_N23-0063_1.sub.fq.gz


# merge the 20% subsample read1 file with the original read2 file
seqtk merge pe sub_reads/HOW_N23-0063_1.sub.fq.gz /data/genomics/workshops/smsc_2024/rawdata/HOW_N23-0063_2.fastq.gz > sub_reads/HOW_N23-0063.interleave.fq


# deinterleave the interleave read file so we get two read files for the subsampled reads: one for read1 and another for read2
sh /data/genomics/workshops/smsc_2024/scripts/deinterleave_fastq.sh < sub_reads/HOW_N23-0063.interleave.fq sub_reads/HOW_N23-0063_1.sub20.fq.gz sub_reads/HOW_N24-0063_2.sub20.fq.gz compress
```

This next script runs MitoZ as a loop over all your samples via their subreads. MitoZ generates a de novo mitogenome assembly and annotates the resulting mitogenome.

```
for sample in 
do
MitoZ.py all --genetic_code 2 --clade 'Chordata' --outprefix ${sample} --thread_number 30 --fastq1 sub_reads/${sample}_1.sub.fq.gz --fastq2 sub_reads/${sample}_2.sub.fq.gz --fastq_quality_shift --fastq_read_length 125 --duplication --insert_size 350 --run_mode 2 --filter_taxa_method 1 --requiring_taxa 'Aves' --species_name 'Hypotaendia owstoni' &> ${sample}.mitoz.log
done
```

### Finding a mitogenone in an assembly using BLASTN

```
cp /data/genomics/workshops/smsc_2024/mitogenomes/*.fa . mitogenomes/
```

```
makeblastdb -in mitogenomes/Hypotaenidia_okinawae.fasta -dbtype nucl
blastn -db Hypotaenidia_okinawae.fasta -query ref_gemome/bHypOws1_hifiasm.bp.p_ctg.fasta -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 12 | sort -V -k12,12r | head -n1 | cut -f1 > GuamRail_mitoscaf.txt
seqtk suubseq ref_genome/bHypOws1_hifiasm.bp.p_ctg.fasta GuamRail_mitoscaf.txt > GuamRail_mitogenome.fasta
```

### Mitogenome Phylogenomics

Now that we have a mitogenome for our samples we are going to do a whole genome phylogeny using mafft and iqtree.

Typically you would want to use the annotations to only include the 12 protein-coding genes (PCGs) but for today we are going to just align the whole mitogenome. Commands on how to concatenate the PCGs can be found: https://github.com/rtfcoimbra/Coimbra-et-al-2021_CurrBiol/blob/main/mitogenomes_workflow.txt

First change into your smsc_2024 directory in the /scratch/genomics/<username>/ and then make a mitogenomes directory and copy the mitogenome fastas

```bash

cd mitogenomes/
cat *.fa > Rail_mitogenomes.fasta

mafft --thread 10 --adjustdirection Rail_mitogenomes.fasta > Rail_mitogenomes.aln
iqtree -nt 10 --pre Rail_mito -s Rail_mitogenomes.aln -B 1000
```

### How to run on an HPC

```bash
#!/bin/sh
#---------------Parameters----------------------#
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=4G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N ANGSD_realSFS_job
#$ -o ANGSD_realSFS_job.log
#---------------Modules-------------------------#
module load bioinformatics/angsd/0.921
module load other_required_modules
#---------------Your Commands-------------------#
echo + date job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME

#Set variables

input_bam_file="path/to/your/input_bam_file.bam"
ancestral_fasta_file="path/to/your/ancestral_fasta_file.fasta"
reference_fasta_file="path/to/your/reference_fasta_file.fasta"
output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

Loop through scaffolds 1 to 19

for i in {1..18}; do
# Run ANGSD command
angsd -P <threads> -i ${input_bam_file} -anc ${ancestral_fasta_file} -dosaf <dosaf_value> -gl <genotype_likelihood_method> -C <base_quality_adjustment> -minQ <min_base_quality> -minmapq <min_mapping_quality> -fold <fold_value> -out ${output_directory}/$SAMPLE.scaffold${i} -ref ${reference_fasta_file} -r HiC_scaffold_${i}

# Run realSFS command
realSFS -nsites <number_of_sites> ${output_directory}/$SAMPLE.scaffold${i}.saf.idx > ${output_directory}/$SAMPLE.scaffold${i}.est.ml

done

echo = date job $JOB_NAME done
```

```bash
#!/bin/sh
#---------------Parameters----------------------#
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=4G,h_data=1G,h_vmem=1G
#$ -cwd
#$ -j y
#$ -N Annotate_EST_ML_job
#$ -o Annotate_EST_ML_job.log
#---------------Modules-------------------------#
#Load any required modules if necessary (uncomment the following line and replace 'module_name')
#module load module_name
#---------------Your Commands-------------------#
echo + date job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#Set variables

output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

#Loop through scaffolds 1 to 19

for i in {1..19}; do
# Add sample name and scaffold number to each line of the output file
awk -v sample="$SAMPLE" -v scaffold="$i" '{print sample, "scaffold" scaffold, $0}' ${output_directory}/$SAMPLE.scaffold${i}.est.ml > ${output_directory}/$SAMPLE.scaffold${i}.est.ml.annotated

# Optional: Move the annotated file to the original file
mv ${output_directory}/$SAMPLE.scaffold${i}.est.ml.annotated ${output_directory}/$SAMPLE.scaffold${i}.est.ml

done

echo = date job $JOB_NAME done
```
