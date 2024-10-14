# PCA and Admixture 

PCA and Admixture are both ways to investigate population structure. 

## Load data

First, we will investigate population structure using principal component analysis. 

- PCA aims to identify the main axes of variation in a dataset with each axis being independent of the next. The first component summarizes the major axis variation and the second the next largest, and so on, until all the cummulative variation is explained. 
- For genomic data, PCA summarizes the major axes of variation in allele frequencies 

To perform PCA, we will use plink, specifically version... 

## Setting us up for analyses

First, lets make a new directory for our analyses. 

```
mkdir /scratch/genomics/YOURNAME/smsc_2024/pop_structure/
mkdir /home/YOURNAME/smsc_2024/day5
```

Next, we will copy all of the already written scripts into this folder so we can slightly alter it: 
```
cp /data/genomics/workshops/smsc_2024/scripts_Day4/*.sh /home/YOURNAME/smsc_2024/day5
```

## Filtering our raw VCF file 

### Check for missingness per individual

First, we need to check out raw VCF file and filter our raw VCF file to contain high quality SNPs. That is, SNPs that we are confident that are correct. 


To do this, we will use VCFtools, which is a very common program used to filter variants and samples in our dataset. 


First,  we will check to see the missingness per individual in our raw VCF file. 

To do this, we will work with the `01_checkmissingness.job` file. 


Submit the script: 
```
qsub 01_checkmissingness.job
```
Let's look at the file to go over what the script is doing: 
```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N missingness                                                                                                                                                              
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
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l.vcf                                                                   
#    

vcftools --vcf ${VCF} \
--missing-indv
echo = `date` job $JOB_NAME done

```

This will output a file called `imiss.out`

Let's open that file to look at the results: 
```
less imiss.out
```
Here's the output: 
```
INDV    N_DATA  N_GENOTYPES_FILTERED    N_MISS  F_MISS                                                                     
Atlantisia_rogersi      10741317        0       2706    0.000251924                                                        
Fulica_atra     10741317        0       175     1.62922e-05                                                                
Grus_americana  10741317        0       9061    0.000843565                                                                
GuamRail_HOW_N23_0063   10741317        0       10      9.30985e-07                                                        
GuamRail_HOW_N23_0068   10741317        0       14      1.30338e-06                                                        
Heliornis_fulica        10741317        0       35441   0.0032995                                                          
Laterallus_jamaicensis_coturniculus     10741317        0       28803   0.00268151                                         
Porphyio_hochstetteri   10741317        0       11989   0.00111616                                                         
Rallus_limicola 10741317        0       2095    0.000195041                                                                
Zapornia_atra   10741317        0       175     1.62922e-05                                                                

```

All of our samples have at least 99% of the some SNPs information at each position of the VCF file. This is great! 

### Filtering to keep best quality SNPs


Next, we will filter the raw VCF to keep the best quality SNPs and remove indels. We will use the script below with VCFtools. 

Open the file:
```
vim 02_FilterSNPs.sh
```
Change YOURNAME to yourname:
```
vim 01_vcftoolsfilter.job
press i
change YOURNAME to your name
Shift x
```

Submit the script: 
```
qsub 02_vcftoolsfilter.job
```

Let's see what the script is doing:
 
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
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l.vcf                                                                      
OUT=/scratch/genomics/YOURNAME/smsc_2024/pop_structure/gatk_allsamples_ptg000001l_filtered                                                                  
#    

vcftools --vcf ${VCF} --minQ 30 --remove-indels --max-missing 0.8 \
--remove-indv Grus_americana --remove-indv Fulica_atra --remove-indv Heliornis_fulica --recode --recode-INFO-all \
--out ${OUT}

```

Specifically, these filters are

- minQ 30 --> keep sites that have at least a minimum quality score of 30. This means we think think the base quality accuracy at this position is 99.9% correct. 
- remove-indels --> this removes indels, which are short insertions and deletions. For hte rest of our analyses, we are just interested in SNPs. 
- max-missing 0.8 --> keep positions where there is data for at least 80% of individuals 


This step will take a couple of minutes. 
Let's compare the output files from the raw VCF and our filtered VCF: 
```
From our unfiltered dataset, we had 10741317 SNPs, as we can see from the imiss file. 

less /data/genomics/workshops/smsc_2024/Results/filterVCF.log 

It states "After filtering, kept 8740710 out of a possible 10741317 Sites"
```
We can see the filtering worked! 



### Use PLINK to remove linked sites

Plink (https://zzz.bwh.harvard.edu/plink/) is a commonly used program to do conduct various population genomic analyses. Mainly, people use it for doing PCA, make input files for other programs, and quantify relatedness. 


An assumption of PCA is that we use independent data -- that is no spurious correlations amont the measured variables. 
- For genomic data, allele frequencies can be correlated due to physical linkage 

First we will open the plink script and change to your name:

```
vim 03_removeLD.job
press i
change YOURNAME to your name
Shift x
```

Then we will submit the job: 
```
qsub 03_removeLD.job
```
This the script below to see what's its doing: 

```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N plinkLD                                                                                                                                                              
#$ -o plinkLD.log                                                                                                                                          
#----------------Modules------------------------- #                                                                                     

module load bio/plink/1.90b7.2
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l_filtered.recode.vcf                                                                    
OUT=/scratch/genomics/YOURNAME/smsc_2024/pop_structure/gatk_allsamples_ptg000001l_filtered_LD

FINAL=/scratch/genomics/YOURNAME/smsc_2024/pop_structure/gatk_allsamples_ptg000001l_filtered_prunedLD
#    

plink --vcf ${VCF} --indep-pairwise 50 5 0.5 --out ${OUT} --const-fid 0 --allow-extra-chr


plink --vcf ${VCF}  --extract ${OUT}.prune.in  --out ${FINAL} --recode --const-fid 0 --dog --allow-extra-chr


echo = `date` job $JOB_NAME done

```

Let's look at the plink output files to understand the file structure: 

```
less /data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD.map

less -S /data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD.ped

```


## PCA!!!!!!! 

Woo! Now we can make our PCA :) We will use plink to infer the PCA. 


First we will open the plink script and change to your name:

Then we will submit the job: 
```
qsub 04_runPCA.job
```
This the script below to see what's its doing: 
```
# /bin/sh                                                                                                                                                                 
# ----------------Parameters---------------------- #                                                                                                                      
#$ -S /bin/sh                                                                                                                                                             
#$ -q sThC.q                                                                                                                                                              
#$ -l mres=10G,h_data=5G,h_vmem=5G                                                                                                                                        
#$ -cwd                                                                                                                                                                   
#$ -j y                                                                                                                                                                   
#$ -N runPCA                                                                                                                                                              
#$ -o runPCA.log                                                                                                                                          
#----------------Modules------------------------- #                                                                                     

module load bio/plink/1.90b7.2
# ----------------Your Commands------------------- #                                                                                                                      
#                                                                                                                                                                         
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                             
#                                                                                                                                                                         
# $SGE_TASK_ID is the value of the task ID for the current task in the array                                                                                              
# write out the job and task ID to a file                                                                                                                                 
#                                                                                                                                                                         
FILE=/data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD_rails                                                                    
#    

plink --file ${FILE} --pca var-wts  --const-fid 0 --allow-extra-chr


echo = `date` job $JOB_NAME done

```

Let's look at the output: 

```
less plink.eigenvec #Principal components
less plink.eigenval #tells us the contribution of each PC in explaining the variation
```

We will now open this in R to plot. Put the file "PCA_rails.csv" into your R directory
```
#in R
install_packages("ggplot2")
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)

#Read in data 
dat <- read.csv ("PCA_rails.csv", header=TRUE)

#plot graph                                      
p <- ggplot(dat,aes(x=PC1,y=PC2,label=Name, color=Name)) + geom_point(size = 3) + theme_classic()  +
  geom_text_repel()
p

#plot Eigenvalues 

Eigenvalues <- c(0.44, .28, .18, .08, .013, .004, -0.0018)
barplot(Eigenvalues)


```

# Individual Admixture Proportions with Admixture 

Now we can infer individual admixture proportions with Admixture. This is a genetic clustering program to define populations and assign individuals to them.  

We will use the plink files we just make to run Admixture. 

```
qsub 05_Admixture.job
```

Here is the script: 

```

# /bin/sh                                                                                                                                                               
# ----------------Parameters---------------------- #                                                                                                                    
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 5                                                                                                                                                        
#$ -q lThM.q                                                                                                                                                            
#$ -l mres=10G,h_data=5G,h_vmem=5G,himem                                                                                                                                
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N Admixtureloop                                                                                                                                                     
#$ -o Admixtureloop.log                                                                                                                                                 
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module load bio/admixture/1.3.0                                                                                                                                         
module load bio/plink/1.90b7.2                                                                                                                                          
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
#                                                                                                                                                                       
INPUT=/data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD_rails                                                                     
OUTPUT=/data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD_rails_bed                                                                
#                                                                                                                                                                       
#plink --file ${INPUT} --make-bed --allow-extra-chr --const-fid 0 --out ${OUTPUT}                                                                                       
#                                                                                                                                                                       
FILE=/data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD_rails_bed.bed                                                              
                                                                                                                                                                        
for i in {2..5}                                                                                                                                                         
do                                                                                                                                                                      
admixture --cv ${FILE} $i > log${i}.out                                                                                                                                 
done                                                                                                                                                                    
#                                                                                                                                                                       
echo = `date` job $JOB_NAME done                                                                                                                                        
#                                                                                                                                                                       

```

This takes a while to finish. We can look at the output files here: 
```
less /data/genomics/workshops/smsc_2024/Results/gatk_allsamples_ptg000001l_filtered_prunedLD_rails_bed.2.Q
```

And now we can plot the data: 

```
library(tidyverse)
library(dplyr)
library (ggplot2)

###################
# Results at K=2 ##
###################
dat <-read.csv ("Admixture_K2.csv", header=TRUE)

## Organize the dataset
data_long <- gather(dat, Admixture, Percentage, Ancestry1:Ancestry2, factor_key=TRUE)

# Plot admixture at K=2
q <- ggplot(data_long, aes(fill=Admixture, y=Percentage, x=Sample)) + 
    geom_bar(position="stack", stat="identity") +scale_fill_manual(values = c("mediumseagreen","maroon4", "orange2", "royalblue4", "red", "honeydew4"))  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
q

###################
# Results at K=3 ##
###################
dat <-read.csv ("Admixture_K3.csv", header=TRUE)

## Organize the dataset
data_long <- gather(dat, Admixture, Percentage, Ancestry1:Ancestry3, factor_key=TRUE)

# Plot admixture at K=2
q <- ggplot(data_long, aes(fill=Admixture, y=Percentage, x=Sample)) + 
    geom_bar(position="stack", stat="identity") +scale_fill_manual(values = c("mediumseagreen","maroon4", "orange2", "royalblue4", "red", "honeydew4"))  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
q

###################
# Results at K=4 ##
###################
dat <-read.csv ("Admixture_K4.csv", header=TRUE)

## Organize the dataset
data_long <- gather(dat, Admixture, Percentage, Ancestry1:Ancestry4, factor_key=TRUE)

# Plot admixture at K=4
q <- ggplot(data_long, aes(fill=Admixture, y=Percentage, x=Sample)) + 
    geom_bar(position="stack", stat="identity") +scale_fill_manual(values = c("mediumseagreen","maroon4", "orange2", "royalblue4", "red", "honeydew4"))  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
q

###################
# Results at K=5 ##
###################
dat <-read.csv ("Admixture_K5.csv", header=TRUE)

## Organize the dataset
data_long <- gather(dat, Admixture, Percentage, Ancestry1:Ancestry5, factor_key=TRUE)

# Plot admixture at K=4
q <- ggplot(data_long, aes(fill=Admixture, y=Percentage, x=Sample)) + 
    geom_bar(position="stack", stat="identity") +scale_fill_manual(values = c("mediumseagreen","maroon4", "orange2", "royalblue4", "red", "honeydew4"))  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
q

```

And we're done!
