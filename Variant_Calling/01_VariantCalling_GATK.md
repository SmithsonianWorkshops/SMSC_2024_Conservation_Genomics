# The Genome Analysis Toolkit (GATK)

The GATK (Genome Analysis Toolkit) is one of the most used programs for genotype calling in sequencing data in model and non model organisms. However, the GATK was designed to analyze human genetic data and all its pipelines are optimized for this purpose. 


[toc]


## Background on file types

- A fastq.gz file contains the raw reads that are not aligned to the reference genome. Let's see what a fastq.gz file looks like: 
```
zless -S /data/genomics/workshops/smsc_2024/rawdata/HOW_N23-0063_1.fastq.gz
```
- A bam file is called a Binary Alignment Map and it is a file where the reads of an individual are aligned to a reference genome. So we know where each read goes on the genome. 
```
module unload gcc/8.5.0
module load bio/samtools 

samtools view /data/genomics/workshops/smsc_2024/BAMS/Fulica_atra.realigned.bam | less -S 
```

- A variant call format (VCF) file contains the variants along the genome for a single individual, or multiple individuals. 
```
less -S /data/genomics/workshops/smsc_2024/VCF/gatk_GuamRails_autosomes_header.vcf
```


## Steps for Alignment to Genotype Calling 

- Map to reference (BAMS) -> HaplotypeCaller -> GenomeDPImport -> Joint call with GenotypeGVCF 

### HaplotypeCaller
HaplotypeCaller calls SNPs and indels via local de-novo assembly of haplotypes. This means when it encounters a region showing signs of variation, it discards mapping information and reassembles the reads in that region. This just allows HaplotypeCaller to be more accurate wen calling regions that are typically difficult to call. 

First, we need to create a new directory: 
```
mkdir /scratch/genomics/YOURNAME/variant_calling/GVCFs
mkdir /home/YOURNAME/day3
```

Copy the template file so we can later it to make our script: 
```
cp /data/genomics/workshops/smsc_2024/jobs/template_file.job /home/hennellyl/day3
```

Then we can run the script:

```     
# /bin/sh                                 
# ----------------Parameters---------------------- #              
#$ -S /bin/sh
#$ -pe mthread 5 
#$ -q lThM.q 
#$ -l mres=10G,h_data=5G,h_vmem=5G,himem  
#$ -cwd 
#$ -j y 
#$ -N haplotypecaller_HOW_N23_0063.job 
#$ -o haplotypecaller_HOW_N23_0063.log                                                                                                                                     
#                                                                                                                                                                         
# ----------------Modules------------------------- #                                                                                                                      
module load bioinformatics/gatk/4.5.0.0                                                                                                                                   
#                                                                                                                                                                         
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
echo + NSLOTS = $NSLOTS                                                                                                                                                   
#                                                                                                                                                                                                                                                               
BAM=/data/genomics/workshops/smsc_2024/BAMS/HOW_N23-0063.realigned.bam                                                                                     
REF=/scratch/genomics/hennellyl/variant_calling/bHypOws1_hifiasm.bp.p_ctg_copy.fasta                                                                                      
#                                                                                                                              
rungatk HaplotypeCaller \                                                                                                                                                 
-I ${BAM} \                                                                                                                                                               
-R ${REF} \                                                                                                                                                               
-ERC GVCF \                                                                                                                                                               
-O /scratch/genomics/YOURNAME/variant_calling/GVCFs/GuamRail_HOW_N23_0063.realigned_duplMarked_SM 

echo = `date` job $JOB_NAME done

```

Let's look at the log file: 
```
less -S haplotypecaller_HOW_N23_0063.log
```

And let's take a look at the GVCF file: 
```
less /data/genomics/workshops/smsc_2024/GVCF/GuamRail_HOW_N23_0063.realigned_duplMarked_SM
```

- Let's look at the header. We have a list of contigs, CHROM, POS, the ref and alternative allele, etc.  


### GenomeDBImport

This step creates a GenomicsDB datastore to speedup joint genotyping. This in prep for the Genotype calling.

First we need to create the directories: 
```
mkdir /scratch/genomics/YOUNAME/variant_calling/GenomeImport_GuamRails 

mkdir /scratch/genomics/hennellyl/variant_calling/GenomeImport_scratch_GuamRails
```
Now we will use an array to run make the GenomeImport folders for each chromosome AT THE SAME TIME. 

- Arrays are super useful for doing work on the cluster more efficiently. 
- What's it doing: for example we have 68 chromosomes in the Guam Rail. We could one a single job file to run an analysis across the whole genome. The array is allowing us to send jobs of the analysis on each chromosome at the same time, thus speeding the analysis by 68x.


/home/hennellyl/jobs

Lets take a look at the files we need for the array: 
```
less  /data/genomics/workshops/smsc_2024/files/autosomes_list.txt
```

We can also take a look at the sample map: 
```
less /data/genomics/workshops/smsc_2024/files/bamlistgaum_samplemap.txt 
```

Here is the script:

```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N GenomicsDBImport                                                                                                                                                    
#$ -o GenomicsDBImport_$JOB_ID-$TASK_ID.log                                                                                                                                       
#$ -t 1-68 -tc 30                                                                                                                                                         
#----------------Modules------------------------- #                                                                                                                       
module load bioinformatics/gatk/4.5.0.0                                                                                                                                   
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
CHR=$(sed -n "${SGE_TASK_ID}p"  /data/genomics/workshops/smsc_2024/files/autosomes_list.txt | cut -f1)                                                                                       
SAMPLEMAP=/data/genomics/workshops/smsc_2024/files/bamlistgaum_samplemap.txt                                                                                                                          
OUTDIR=/scratch/genomics/YOURNAME/variant_calling/GenomeImport_GuamRails                                                                                                 
TEMPDIR=/scratch/genomics/YOURNAME/variant_calling/GenomeImport_scratch_GuamRails                                                                                        
#                                                                                                                                                                         
rungatk GenomicsDBImport \                                                                                                                                                
--genomicsdb-workspace-path ${OUTDIR}/${CHR}_gvcf_db \                                                                                                                    
--batch-size 50 \                                                                                                                                                         
-L ${CHR} \                                                                                                                                                               
--sample-name-map ${SAMPLEMAP} \                                                                                                                                          
--tmp-dir ${TEMPDIR}                                                                                                                                                      
                                                                                                                                                                          
                                                                                                                                                                          
echo = `date` job $JOB_NAME done                                                                                                                                          

```

Let's see how its working. We can cd into the GenomeDBImport: 

```
cd /scratch/genomics/YOURNAME/variant_calling/GenomeImport_GuamRails
```

We should see the GenomeBPImport start making the databases by chromosome
### GenotypeGVCFs

Now we can perform joint genotyping that were pre-called with HaplotypeCaller. 
- It will look at the available information for each site from both variant and non-variant alleles across all samples, and then produce a multi-sample VCF that contains only sites that are found as variant in at least one sample. 

https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs

Let's make a directory for our VCF output: 
```
mkdir /scratch/genomics/YOURNAME/variant_calling/GenotypeGATK_bothGuamRail

mkdir /scratch/genomics/YOURNAME/variant_calling/GenotypeGATK_scratch_bothGuamRail/
```
Now let's run the job file:
```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N GenotypeGATK                                                                                                                                                        
#$ -o genotypegatk_$JOB_ID-$TASK_ID.log                                                                                                                                  
#$ -t 1-68 -tc 30                                                                                                                                                         
#----------------Modules------------------------- #                                                                                                                       
module load bioinformatics/gatk/4.5.0.0                                                                                                                                   
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
CHR=$(sed -n "${SGE_TASK_ID}p"  /data/genomics/workshops/smsc_2024/files/autosomes_list.txt | cut -f1)                                                                                       
REF=/scratch/genomics/hennellyl/variant_calling/bHypOws1_hifiasm.bp.p_ctg_copy.fasta                                                                                      
OUTDIR=/scratch/genomics/YOURNAME/variant_calling/GenotypeGATK_bothGuamRail                                                                                              
TEMPDIR=/scratch/genomics/YOURNAME/variant_calling/GenotypeGATK_scratch_bothGuamRail/                                                                                    
#                                                                                                                                                                         
cd /scratch/genomics/YOURNAME/variant_calling/GenomeImport_GuamRails                                                                                                               
#                                                                                                                                                                         
rungatk GenotypeGVCFs \                                                                                                                                                   
-R ${REF} \                                                                                                                                                               
-V gendb://${CHR}_gvcf_db \                                                                                                                                               
-O ${OUTDIR}/gatk_GuamRails_${CHR}.vcf.gz \                                                                                                                           
--tmp-dir ${TEMPDIR}                                                                                                                                                      
#                                                                                                                                                                         
#                                                                                                                                                                         
echo = `date` job $JOB_NAME done                               
```


Let's look at the multi-sample VCF:
```
less -S /data/genomics/workshops/smsc_2024/VCF/gatk_GuamRails_autosomes_header.vcf
```
- We can see the FORMAT column, which as ...
- Under the two bird individuals, we can the information associated with each genotype call. For example: 
```
0|1:4,2:6:61:0|1:16261_A_C:61,0,162:16261
```
- The `0|1` indicates the 0 is the reference allele, the 1 is the alternative allele, the vertical pipe | is showing a phased genotype, and the slash / indicates the chromosome of the alleles is unknown. 
- The `4,2` indicates the allele count at that site for the REF and ALT allele. Here, we had 6 reads that overlapped at this genomic position. Among those 6 reads, there were 4 that were Reference allele, and 2 that were Alternative allele. Thus, the program called it as a heterozygote. 



### VariantFiltration

Now we will filter our variants based on certain criteria. 

https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration

Meaning of filters below: 
- filter "DP < 1800" -> remove sites that have very high depth (why is this important?)
- filter FS -> FisherStrand, accounts for strand bias at the site. Strand bias is whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. This results in the incorrect amount of evidence 



```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N Filter                                                                                                                                                              
#$ -o variant_filter_ID.log                                                                                                                                          
#----------------Modules------------------------- #                                                                                                                       
module load bioinformatics/gatk/4.5.0.0                                                                                                                                   
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_GuamRails_autosomes_header.vcf

OUT=/scratch/genomics/YOURNAME/variant_calling/GenotypeGATK_bothGuamRail/gatk_GuamRails_autosomes_filtered.vcf                                                                    
#                                                                                                                                                                         
#                                                                                                                                                                         
rungatk VariantFiltration \                                                                                                                                               
-V ${VCF} \                                                                                                                                                               
-filter "DP < 1800" --filter-name "DP1800" \                                                                                                                                
-filter "SOR > 3.0" --filter-name "SOR3" \                                                                                                                                
-filter "FS > 60.0" --filter-name "FS60" \                                                                                                                                
-filter "MQ < 40.0" --filter-name "MQ40" \                                                                                                                                
-filter "MQRankSum < -12.5"  --filter-name "MQRankSum-12.5" \                                                                                                             
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \                                                                                                        
-O ${OUT}                                                                                                                                                                 
#                                                                                                                                                                         
#                                                                                                                                                                         
echo = `date` job $JOB_NAME done    
```

### VCFtools Filtering
Introduction to VCFtools (https://vcftools.sourceforge.net/man_latest.html)

We will use this program a lot in population genomics. It is common for filtering VCF datasets. 


The script: 
```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N Filter                                                                                                                                                              
#$ -o vcftoolsfilter.log                                                                                                                                          
#----------------Modules------------------------- #                                                                                                                       
module unload gcc/8.5.0
module load bio/vcftools
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_GuamRails_autosomes_header.vcf                                                                      
OUT=/scratch/genomics/YOURNAME/variant_calling/GenotypeGATK/gatk_Guam_N23_0063_autosomes_filtered                                                                  
#    

vcftools --vcf ${VCF} \
--minQ 30 --remove-indels \
--recode --recode-INFO-all \
--out ${OUT}
```
Let's look at the log file:
```
less vcftoolsfilter.log
```

To filter our VCF dataset by: 
- no indels (--remove-indels)
    - what are indels? It stands for Insertion (In) or Deletion (Del). These are structural variants, which means there is an insertion in the genome, which alters the structure of the genome.
- keep sites that have a minimum quality score of 30
    - The genotype quality score represents the confidence that a genotype assignment is correct. A low GQ indicates there isn't enough evidence to confidently choose one genotype over another. 
        - A GQ of 30 indicates virtually all reads have no errors. A GQ of 20 corresponds to an error rate of 1%. GQ of 10 indicates a estimated error rate of 10%. 

- Mapping Quality Score (minQ)- quantify the probability that a read is misplaced. An example is many repeated areas of the genome have low mapping quality score. 

## Check the missingness

We can check the missingness per individual. A sample with high missingness means that it has no genotypes present in the VCF file at a large percentage of the genome. 

```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N Filter                                                                                                                                                              
#$ -o missingness.log                                                                                                                                          
#----------------Modules------------------------- #                                                                                                                       
module unload gcc/8.5.0
module load bio/vcftools
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                        
VCF=/scratch/genomics/YOURNAME/variant_calling/GenotypeGATK/gatk_Guam_N23_0063_autosomes_filtered.recode.vcf                                                                   
#    

vcftools --vcf ${VCF} \
--missing-indv

```

Now we can look at the output:
```
less out.miss 
```



